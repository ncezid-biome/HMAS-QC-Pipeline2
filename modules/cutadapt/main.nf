process cutadapt {
    publishDir "${params.final_outdir}/${sample}", mode: 'copy', pattern: "*.fastq"
    tag "${sample}"
    cpus = "${params.maxcpus}"
    memory = "${params.medmems}"
    errorStrategy 'retry'
    maxRetries 3
    // debug true

    maxForks = "${params.maxcutadapts}"

    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path ("cutadapt/${sample}.trimmed.1.fastq"), path ("cutadapt/${sample}.trimmed.2.fastq"), emit: cutadapt_fastq, optional: true
    tuple val(sample), path ("cutadapt/${sample}.cutadapt.json"), emit: cutadapt_json, optional: true

    shell:
    '''
    mkdir -p cutadapt
    run_cutadapt.py -f !{reads[0]} -r !{reads[1]} \
                    -o cutadapt -s !{sample} -p !{params.primer} \
                    -e !{params.cutadapt_maxerror} -l !{params.cutadapt_minlength} \
                    -t !{params.cutadapt_thread} -b !{params.cutadapt_long}

    '''

}