process make_count_table {
    publishDir "${params.final_outdir}/${sample}/temp", pattern: "*.count_table", mode: 'copy'
    publishDir "${params.final_outdir}/${sample}", pattern: "*.csv", mode: 'copy'
    tag "${sample}"
    // debug true
    cpus = "${params.mincpus}"
    memory = "${params.maxmems}"

    input:
    tuple val(sample), path (match_file), path (fasta_file), path (ch_primer_file)

    output:
    path ("${sample}.final.count_table"), emit:count, optional:true
    path ("${sample}.csv"), emit:report, optional:true
    path ("${sample}.primer_stats.tsv"), emit:primer_stats, optional:true
    path ("${sample}.read_length.tsv"), emit:read_length, optional:true

    shell:
    '''
    if [ -s !{match_file} ]; then
        make_count_table.py -o !{sample}.final.count_table -m !{match_file}
        create_report.py -s !{sample} -c !{sample}.final.count_table -p !{ch_primer_file} -o !{sample}.csv \
                         -q !{sample}.primer_stats.tsv -f !{fasta_file} -l !{sample}.read_length.tsv
    else
        echo "!{match_file} is empty !"
    fi
    '''

}