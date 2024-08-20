process pair_merging {
    publishDir "${params.final_outdir}/${sample}/temp", mode: 'copy', pattern: "*.fastq"
    tag "${sample}"
    cpus = "${params.medcpus}"
    errorStrategy 'ignore'

    input:
    tuple val(sample), path(reads1), path(reads2)

    output:
    tuple val(sample), path ("${sample}.fastq"), emit: fastq, optional: true
    tuple val(sample), path ("${sample}_pear.log"), emit: log, optional: true
    path ("${sample}_pear_log.csv"), emit: log_csv, optional: true

    shell:
    '''
    if [ -s !{reads1} ] && [ -s !{reads2} ]; then
        pear -f !{reads1} -r !{reads2} -o !{sample} -q !{params.merging_minquality} \
                                                -m !{params.merging_maxlength} \
                                                -v !{params.merging_minoverlap} -j !{{params.medcpus}} \
                                                > !{sample}_pear.log
        mv !{sample}.assembled.fastq !{sample}.fastq
        parse_pear_log.py -p !{sample}_pear.log -s !{sample} -o !{sample}_pear_log.csv
    else
        echo "either !{reads1} or !{reads2} is empty !"
    fi

    '''
}