process split_by_adapter {
    publishDir "${params.final_outdir}/${sample}/by_adapter", mode: 'copy', pattern: "by_adapter/*"
    tag "${sample}"
    cpus = "${params.mincpus}"
    memory = "${params.medmems}"
    errorStrategy 'retry'
    maxRetries 3

    maxForks = "${params.maxcutadapts}"

    when:
    params.split_by_adapter

    input:
    tuple val(sample), path(matched_r1), path(matched_r2)

    output:
    path("by_adapter/*.fastq"), optional: true

    shell:
    '''
    mkdir -p by_adapter
    bin_fastq_by_adapter.py \
        -r1 !{matched_r1} \
        -r2 !{matched_r2} \
        -o by_adapter
    '''
}
