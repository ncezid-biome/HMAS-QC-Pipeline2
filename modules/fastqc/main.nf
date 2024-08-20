process FASTQC {
    publishDir "${params.final_outdir}/${sample}/fastqc", mode: 'copy'
    tag( "${sample}" )

    input:
    tuple val(sample), path(reads)
 
    output:
    tuple val(sample), path ("*_fastqc.{zip,html}"), emit: fastqc_results
 
    script:
    """
    fastqc -q $reads
    """
}