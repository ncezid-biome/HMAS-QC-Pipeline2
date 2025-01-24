#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define the pipeline version
def pipeline_version = '1.2.1' // Replace with actual version or load dynamically

// Define the timestamp
def timestamp = new Date().format("yyyyMMdd_HHmmss")

// Define the output directory with the version and timestamp at runtime
params.final_outdir = params.outdir ? "${params.outdir.replaceAll('/+$', '')}_v${pipeline_version}_${timestamp}" : "hmas2_results_v${pipeline_version}_${timestamp}"
params.file_extension = "_v${pipeline_version}_${timestamp}"

Channel
    // search for pair-end raw reads files in the given folder or any subfolders
   .fromFilePairs(["${params.reads}/*_R{1,2}*.fastq.gz", "${params.reads}/**/*_R{1,2}*.fastq.gz"], size: 2)
//   .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
    .map{ reads -> tuple(reads[0].replaceAll(~/_L[0-9]+/,""), reads[1]) }
  .set { paired_reads }

Channel.fromPath(params.multiqc_config, checkIfExists: true).set { ch_config_for_multiqc }
Channel.fromPath(params.custom_logo, checkIfExists: true).set { ch_logo_for_multiqc }
Channel.fromPath(params.primer, checkIfExists: true).set { ch_primer_file }


include { FASTQC as FASTQC_RAW } from './modules/fastqc/main.nf' 
include { cutadapt } from './modules/cutadapt/main.nf' 
include { pair_merging } from './modules/pair_merging/main.nf' 
include { quality_filtering; dereplication; denoising; search_exact } from './modules/vsearch/main.nf'
// include { hashing } from './modules/local/hash' 
include { combine_reports } from './modules/local/combine_reports.nf'
include { combine_logs as combine_logs_pear } from './modules/local/combine_logs.nf'
include { combine_logs as combine_logs_qfilter } from './modules/local/combine_logs.nf' 
include { combine_logs as combine_logs_derep } from './modules/local/combine_logs.nf' 
include { combine_logs as combine_logs_denoise } from './modules/local/combine_logs.nf' 
include { make_count_table } from './modules/local/make_count_table.nf' 
include { multiqc } from './modules/multiqc/main.nf' 

workflow {
    // Filter out file pairs containing "Undetermined"
    paired_reads = paired_reads.filter { pair -> 
    !new File(pair[0]).getName().toLowerCase().startsWith("undetermined")}

    FASTQC_RAW(paired_reads)
 // removed_primer_reads_ch = cutadapt(paired_reads)
    paired_reads.combine(ch_primer_file).set{ ch_for_cutadapt }
    removed_primer_reads_ch = cutadapt(ch_for_cutadapt)
    merged_reads_ch = pair_merging(removed_primer_reads_ch.cutadapt_fastq)
    filered_reads_ch = quality_filtering(merged_reads_ch.fastq)
    unique_reads_ch = dereplication(filered_reads_ch.fasta)
    denoisded_reads_ch = denoising(unique_reads_ch.fasta)
    // hashing(denoisded_reads_ch.unique)
    before_search_ch = filered_reads_ch.fasta.join(denoisded_reads_ch.unique)
    match_file_ch = search_exact(before_search_ch)
    // collectFile will instead concatenate all the file contents and write it into a single file
    // which is not what we want.  We want to read each file separately, for all the files
    before_count_table_ch = match_file_ch.join(denoisded_reads_ch.unique)
    before_count_table_ch.combine(ch_primer_file).set{ ch_for_make_count_table }
    reports_file_ch = make_count_table(ch_for_make_count_table)
    // reports_file_ch = make_count_table(before_count_table_ch)
    combined_report_ch = combine_reports(reports_file_ch.report.collect(), \
                                         reports_file_ch.primer_stats.collect(), \
                                         reports_file_ch.read_length.collect(), \
                                         ch_primer_file)

    pear_log_ch = combine_logs_pear(merged_reads_ch.log_csv.collect(), Channel.value('pear'))
    qfilter_log_ch = combine_logs_qfilter(filered_reads_ch.log_csv.collect(), Channel.value('qfilter'))
    derep_log_ch = combine_logs_derep(unique_reads_ch.log_csv.collect(), Channel.value('dereplication'))
    denoise_log_ch = combine_logs_denoise(denoisded_reads_ch.log_csv.collect(), Channel.value('denoise'))

    // add fastqc and cutadapt log files (these are existing modules in MultiQC)
    Channel.empty()
        .mix( FASTQC_RAW.out.fastqc_results )
        .mix( removed_primer_reads_ch.cutadapt_json )
        .map { sample, files -> files }
        .collect().ifEmpty([])
        .set { log_files }

    // add custom content log files
    multiqc(log_files
        .combine(ch_logo_for_multiqc)
        .combine(pear_log_ch.log)
        .combine(qfilter_log_ch.log)
        .combine(derep_log_ch)
        .combine(denoise_log_ch)
        .combine(combined_report_ch.primer_stats_mqc)
        .combine(combined_report_ch.read_length_mqc)
        .combine(combined_report_ch.report_mqc), ch_config_for_multiqc)
}
