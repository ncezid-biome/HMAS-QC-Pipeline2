#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.primers = workflow.launchDir + '/M3235_22_024.primers.test'
params.outdir = workflow.launchDir + '/output_samplebase'
// params.reads = workflow.launchDir
params.reads = workflow.launchDir + '/data'
params.oligo = workflow.launchDir + '/data' + '/M3235_22_024.oligos'

Channel
  .fromFilePairs("${params.reads}/*_R{1,2}*.fastq.gz",size: 2)
 .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
//   .view()
  .set { paired_reads }

// Channel
//   .fromFilePairs("${params.reads}/*.{1,2}.fastq",size: 2)
// //   .map{ reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
// //    .view()
//   .set { paired_reads2 }

// Channel
//     .fromPath(params.primers)
//     .splitCsv(header:false, sep:"\t")
//     .set{primer_ch}


process cutadapt {
    // publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus 20
    memory = 2.GB
    // container 'dceoy/cutadapt:latest'
    // debug true
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path ("cutadapt/${sample}*.1.fastq"), path ("cutadapt/${sample}*.2.fastq")

    shell:
    '''
    mkdir -p cutadapt
    #python !{workflow.projectDir}/bin/run_cutadapt.py -f !{reads[0]} -r !{reads[1]} \
    #                                 -o cutadapt -s !{sample} -p !{params.oligo}
    run_cutadapt.py -f !{reads[0]} -r !{reads[1]} \
                                     -o cutadapt -s !{sample} -p !{params.oligo}
    '''

}

process concat_reads {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    // debug true

    input:
    tuple val(sample), path (reads1), path (reads2)

    output:
    tuple val(sample), path ("${sample}.1.fastq"), path ("${sample}.2.fastq")

    shell:
    '''
    cat !{reads1} >> !{sample}.1.fastq
    cat !{reads2} >> !{sample}.2.fastq

    '''

}

process pair_merging {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    cpus = 18

    input:
    tuple val(sample), path(reads1), path(reads2)

    output:
    // path ("${sample}.fastq")
    tuple val(sample), path ("${sample}.fastq")

    shell:
    '''
    pear -f !{reads1} -r !{reads2} -o !{sample} -q 26 -m 325 -v 20 -j 20
    mv !{sample}.assembled.fastq !{sample}.fastq

    '''
}

process quality_filtering {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    // tag "quality_filtering"
    tag "${sample}"

    input:
    // path ("*.fastq")
    tuple val(sample), path (fastq)

    output:
    // path ("output.fasta")
    tuple val(sample), path ("${sample}.fasta")

    shell:
    '''
    #cat *.fastq >> output.fastq
    vsearch --fastx_filter !{fastq} --fastq_maxee 1 --fastaout !{sample}.fasta

    #remove space between seq_id and =adapter 
    #choose not to use sed, because it could be very slow
    #sed -r -i 's/\s+adapter=/=/g' output.fasta
    #python !{workflow.projectDir}/bin/remove_space.py -f output.fasta

    remove_space.py -f !{sample}.fasta
    '''

}

process dereplication {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    // tag "dereplication"
    tag "${sample}"

    input:
    tuple val(sample), path (fasta)

    output:
    // path ("unique.fasta")
    tuple val(sample), path ("${sample}.unique.fasta")

    shell:
    '''
    vsearch --derep_fulllength !{fasta} --output !{sample}.unique.fasta --sizeout --relabel_keep

    '''

}

process denoising {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    // tag "denoising"
    tag "${sample}"
    // debug true
    cpus = 3

    input:
    tuple val(sample), path (fasta)

    output:
    tuple val(sample), path ("${sample}.final.unique.fasta")

    shell:
    '''
    vsearch --cluster_unoise !{fasta} --minsize 10 --unoise_alpha 2 --centroids !{sample}.unique.unoise.fasta
    vsearch --uchime3_denovo !{sample}.unique.unoise.fasta --nonchimeras !{sample}.final.unique.fasta

    '''

}


process search_exact {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    // tag "searching for exact seqs"
    tag "${sample}"
    cpus = 9

    input:
    // path (final_unique_fasta)
    // path (output_fasta)
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
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    // tag "generating abundance table"
    tag "${sample}"
    // debug true
    cpus = 4
    memory = 16.GB

    input:
    tuple val(sample), path (match_file)

    output:
    path ("${sample}.final.count_table.*"), optional:true
    path ("${sample}.final.count_table"), optional:true

    shell:
    '''
    #python !{workflow.projectDir}/bin/make_count_table.py -o final.count_table -m !{match_file}
    if [ -s !{match_file} ]; then
        make_count_table.py -o !{sample}.final.count_table -m !{match_file}
    else
        echo "!{match_file} is empty !"
    fi
    '''

}

workflow {
    removed_primer_reads_ch = cutadapt(paired_reads)
    clean_reads_ch = concat_reads(removed_primer_reads_ch)
    // merged_reads_ch = pair_merging(clean_reads_ch).collect()
    merged_reads_ch = pair_merging(clean_reads_ch)
    filered_reads_ch = quality_filtering(merged_reads_ch)
    unique_reads_ch = dereplication(filered_reads_ch)
    denoisded_reads_ch = denoising(unique_reads_ch)

    before_search_ch = filered_reads_ch.join(denoisded_reads_ch)
    // match_file_ch = search_exact(denoisded_reads_ch,filered_reads_ch)
    match_file_ch = search_exact(before_search_ch)
    make_count_table(match_file_ch)
}