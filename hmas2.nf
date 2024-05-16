#!/usr/bin/env nextflow
nextflow.enable.dsl=2


Channel
    // search for pair-end raw reads files in the given folder or any subfolders
   .fromFilePairs(["${params.reads}/*_R{1,2}*.fastq.gz", "${params.reads}/**/*_R{1,2}*.fastq.gz"], size: 2)
 .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .set { paired_reads }

process cutadapt {
    // publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    cpus = "${params.maxcpus}"
    memory = "${params.medmems}"
    errorStrategy 'retry'
    maxRetries 3
    debug true

    maxForks = "${params.maxcutadapts}"

    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path ("cutadapt/${sample}.1.fastq"), path ("cutadapt/${sample}.2.fastq")

    shell:
    '''
    mkdir -p cutadapt
    run_cutadapt.py -f !{reads[0]} -r !{reads[1]} \
                    -o cutadapt -s !{sample} -p !{params.primer} \
                    -e !{params.cutadapt_maxerror} -l !{params.cutadapt_minlength} \
                    -t !{params.cutadapt_thread} -b !{params.cutadapt_long}

    '''

}

process pair_merging {
    publishDir "${params.outdir}/${sample}/temp", mode: 'copy'
    tag "${sample}"
    cpus = "${params.medcpus}"
    errorStrategy 'ignore'

    input:
    tuple val(sample), path(reads1), path(reads2)

    output:
    tuple val(sample), path ("${sample}.fastq"), optional:true

    shell:
    '''
    if [ -s !{reads1} ] && [ -s !{reads2} ]; then
        pear -f !{reads1} -r !{reads2} -o !{sample} -q !{params.merging_minquality} \
                                                -m !{params.merging_maxlength} \
                                                -v !{params.merging_minoverlap} -j !{{params.medcpus}}
        mv !{sample}.assembled.fastq !{sample}.fastq
    else
        echo "either !{reads1} or !{reads2} is empty !"
    fi

    '''
}

process quality_filtering {
    publishDir "${params.outdir}/${sample}/temp", mode: 'copy'
    tag "${sample}"

    input:
    tuple val(sample), path (fastq)

    output:
    tuple val(sample), path ("${sample}.fasta"), optional:true

    shell:
    '''
    vsearch --fastx_filter !{fastq} --fastq_maxee 1 --fastaout !{sample}.fasta

    #remove space between seq_id and =adapter 
    #choose not to use linux sed , because it could be very slow
    remove_space.py -f !{sample}.fasta
    '''

}

process dereplication {
    publishDir "${params.outdir}/${sample}/temp", mode: 'copy'
    tag "${sample}"

    input:
    tuple val(sample), path (fasta)

    output:
    tuple val(sample), path ("${sample}.unique.fasta")

    shell:
    '''
    vsearch --derep_fulllength !{fasta} --output !{sample}.unique.fasta --sizeout --relabel_keep

    '''

}

process denoising {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    // debug true
    cpus = "${params.mincpus}"

    input:
    tuple val(sample), path (fasta)

    output:
    tuple val(sample), path ("${sample}.final.unique.fasta"), optional:true

    shell:
    '''
    vsearch --cluster_unoise !{fasta} --minsize !{params.denoising_minsize} \
                                      --unoise_alpha !{params.denoising_alpha} \
                                      --centroids !{sample}.unique.unoise.fasta
    
    #remove potential chimeras
    vsearch --uchime3_denovo !{sample}.unique.unoise.fasta --nonchimeras !{sample}.final.unique.fasta

    '''

}


process search_exact {
    publishDir "${params.outdir}/${sample}/temp", mode: 'copy'
    tag "${sample}"
    cpus = "${params.medcpus}"

    input:
    tuple val(sample), path (output_fasta), path (final_unique_fasta)

    output:
    tuple val(sample), path ("${sample}.output.match.final.txt")

    shell:
    '''
    vsearch --search_exact !{output_fasta} -db !{final_unique_fasta} \
            --userfields target+query --userout  !{sample}.output.match.final.txt

    '''

}

process make_count_table {
    publishDir "${params.outdir}/${sample}/temp", pattern: "*.count_table", mode: 'copy'
    publishDir "${params.outdir}/${sample}", pattern: "*.csv", mode: 'copy'
    tag "${sample}"
    // debug true
    cpus = "${params.mincpus}"
    memory = "${params.maxmems}"

    input:
    tuple val(sample), path (match_file)

    output:
    path ("${sample}.final.count_table"), emit:count, optional:true
    path ("${sample}.csv"), emit:report, optional:true

    shell:
    '''
    if [ -s !{match_file} ]; then
        make_count_table.py -o !{sample}.final.count_table -m !{match_file}
        create_report.py -s !{sample} -c !{sample}.final.count_table -p !{params.primer} -o !{sample}
    else
        echo "!{match_file} is empty !"
    fi
    '''

}

process combine_reports {
    publishDir "${params.outdir}", mode: 'copy'
    tag "combine reports"
    debug true

    input:
    path (reports_file)

    output:
    path ("report.csv"), optional:true

    shell:
    '''
    #!{reports_file} is passed in as a string (sapce delimited) concatenation of all sample.csv file
    # Ex.  sample1.csv sample2.csv sample3.csv 
    # which will then be split by the script to read each csv file
    #combine_reports.py -o report.csv -p "!{reports_file}"
    combine_reports.py -o report.csv -p "!{reports_file}" -i "!{params.reads}"

    '''

}

workflow {
    // Filter out file pairs containing "Undetermined"
    paired_reads = paired_reads.filter { pair -> 
    !new File(pair[0]).getName().toLowerCase().startsWith("undetermined")}

    removed_primer_reads_ch = cutadapt(paired_reads)
    merged_reads_ch = pair_merging(removed_primer_reads_ch)
    filered_reads_ch = quality_filtering(merged_reads_ch)
    unique_reads_ch = dereplication(filered_reads_ch)
    denoisded_reads_ch = denoising(unique_reads_ch)
    before_search_ch = filered_reads_ch.join(denoisded_reads_ch)
    match_file_ch = search_exact(before_search_ch)
    // collectFile will instead concatenate all the file contents and write it into a single file
    // which is not what we want.  We want to read each file separately, for all the files
    reports_file_ch = make_count_table(match_file_ch).report.collect()
    combine_reports(reports_file_ch)
}

