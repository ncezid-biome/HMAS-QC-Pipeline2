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