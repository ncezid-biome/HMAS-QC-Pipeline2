process hashing {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    // debug true
    cpus = "${params.mincpus}"

    input:
    tuple val(sample), path (fasta)

    output:
    path ("${sample}.final.unique.hashed"), emit:hash, optional:true

    shell:
    '''
    hash_sequence.py -i !{fasta} -o !{sample}.final.unique.hashed

    '''

}