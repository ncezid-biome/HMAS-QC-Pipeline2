#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// LOGGING cofig
def logMessage(msg) {
    def LOG_FILE = "step_mothur_pipeline.log"
    new File("${LOG_FILE}").withWriterAppend { writer ->
        writer.println("[${new Date()}] $msg")
    }
}


// Define the timestamp
def timestamp = new Date().format("yyyyMMdd_HHmmss")

// Define the output directory with the version and timestamp at runtime
params.final_outdir = params.outdir ? "${params.outdir.replaceAll('/+$', '')}_v${params.pipeline_version}_${timestamp}" : "hmas2_results_v${params.pipeline_version}_${timestamp}"
params.file_extension = "_v${params.pipeline_version}_${timestamp}"

// Track occurrences of read names
// updated to guard against name collision in FASTQ files. In case of dupplicate name
// FASTQ files, they'll be appended with _2, _3 etc. (updated on disk as well)
def name_counts = [:]
Channel
    // search for pair-end raw reads files in the given folder or any subfolders
    .fromFilePairs(
        ["${params.reads}/*_R{1,2}*.fastq.gz", "${params.reads}/**/*_R{1,2}*.fastq.gz"],
        size: 2
    )
    .map { reads_name, reads_paths ->
        // Remove _L### from the sample name
        def cleaned_name = reads_name.replaceAll(/_L[0-9]+/, '')

        // If this name has been seen before, increment counter and append suffix
        def count = name_counts.get(cleaned_name, 0) + 1
        name_counts[cleaned_name] = count

        // Append suffix only if it's a duplicate (i.e., count > 1)
        def final_name = (count > 1) ? "${cleaned_name}_${count}" : cleaned_name

        // List to hold the updated paths
        def updated_paths = []

        if (count > 1) {
            reads_paths.each { path ->
                def base_name = path.getName()
                def suffix = count
                def new_name = base_name.replaceAll(cleaned_name, "${cleaned_name}_${suffix}")
                def new_path = path.getParent().resolve(new_name)

                // Keep incrementing suffix if file already exists
                while (new_path.exists()) {
                    suffix += 1
                    new_name = base_name.replaceAll(cleaned_name, "${cleaned_name}_${suffix}")
                    new_path = path.getParent().resolve(new_name)
                }

                // Update name_counts for this new suffix
                name_counts[cleaned_name] = suffix

                // Rename the file
                path.renameTo(new_path)
                println "Renamed: ${path} to ${new_path}"
                updated_paths << new_path
            }
        } else {
            updated_paths = reads_paths
        }

        // Return the final name and the updated paths (either renamed or original)
        tuple(final_name, updated_paths)
    }
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
    logMessage("step_mothur started")
    logMessage("processing reads from: ${params.reads}")
    logMessage("will save results into: ${params.final_outdir}")
    
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

    process make_command_yaml {

        output:
        path "cli_mqc.txt" , optional:true, emit: CLI

        script:
        """
        parse_commandline.py \
            --cmdline "${workflow.commandLine}" \
            --params_str "--reads ${params.reads} --primer ${params.primer}" \
            --output cli_mqc.txt
        """

    }
    make_command_yaml_ch = make_command_yaml()

    // mix all the version files into one channel
    all_versions = Channel
        .empty()
        .mix(merged_reads_ch.versions)
        .mix(filered_reads_ch.versions)
        .mix(combined_report_ch.versions)

    collected_versions = all_versions
        .collectFile(name: 'software_versions.yml')
        .ifEmpty([]) // returns an empty list if nothing is collected


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
        .combine(make_command_yaml_ch.CLI)
        .combine(combined_report_ch.report_mqc), ch_config_for_multiqc, collected_versions)
}

// Capture Nextflow pipeline completion stats and append to the log file
workflow.onComplete {
    logMessage("step_mothur finished!")
    def stats = workflow.stats
    logMessage("  - Processes executed: ${stats.succeedCount}")
    logMessage("  - Processes failed: ${stats.failedCount}")
    logMessage("  - Processes cached: ${stats.cachedCount}")
    logMessage("  - Workflow duration: ${workflow.duration}")
}

// Capture pipeline errors
workflow.onError { 
    logMessage("${workflow.errorReport}")
}
