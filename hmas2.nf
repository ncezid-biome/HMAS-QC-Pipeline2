#!/usr/bin/env nextflow
nextflow.enable.dsl=2


Channel
    // search for pair-end raw reads files in the given folder or any subfolders
   .fromFilePairs(["${params.reads}/*_R{1,2}*.fastq.gz", "${params.reads}/**/*_R{1,2}*.fastq.gz"], size: 2)
 .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .set { paired_reads }


include { FASTQC as FASTQC_RAW } from './modules/fastqc'
include { cutadapt } from './modules/cutadapt'
include { pair_merging } from './modules/pair_merging'
include { quality_filtering; dereplication; denoising; search_exact } from './modules/vsearch'
// include { hashing } from './modules/local/hash'
include { combine_reports } from './modules/local/combine_reports'
include { make_count_table } from './modules/local/make_count_table'
include { multiqc } from './modules/multiqc'

workflow {
    // Filter out file pairs containing "Undetermined"
    paired_reads = paired_reads.filter { pair -> 
    !new File(pair[0]).getName().toLowerCase().startsWith("undetermined")}

    FASTQC_RAW(paired_reads)
    removed_primer_reads_ch = cutadapt(paired_reads)
    merged_reads_ch = pair_merging(removed_primer_reads_ch.cutadapt_fastq)
    filered_reads_ch = quality_filtering(merged_reads_ch)
    unique_reads_ch = dereplication(filered_reads_ch)
    denoisded_reads_ch = denoising(unique_reads_ch)
    // hashing(denoisded_reads_ch.unique)
    before_search_ch = filered_reads_ch.join(denoisded_reads_ch.unique)
    match_file_ch = search_exact(before_search_ch)
    // collectFile will instead concatenate all the file contents and write it into a single file
    // which is not what we want.  We want to read each file separately, for all the files
    reports_file_ch = make_count_table(match_file_ch).report.collect()
    combine_reports(reports_file_ch)

    Channel.empty()
        .mix( FASTQC_RAW.out.fastqc_results )
        .mix( removed_primer_reads_ch.cutadapt_json )
        .map { sample, files -> files }
        .collect().ifEmpty([])
        .set { log_files }

    multiqc(log_files)
}
