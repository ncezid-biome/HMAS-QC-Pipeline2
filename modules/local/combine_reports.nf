process combine_reports {
    publishDir "${params.final_outdir}", mode: 'copy'
    tag "combine reports"
    // debug true

    input:
    path (reports_file)
    path (primer_stats)
    path (read_length)

    output:
    path ("report*.csv"), emit: csv, optional:true
    path ("report_mqc.yaml"), emit: report_mqc, optional:true
    path ("primer_stats_mqc.yaml"), emit: primer_stats_mqc, optional:true
    path ("read_length_mqc.yaml"), emit: read_length_mqc, optional:true

    shell:
    '''
    #!{reports_file} is passed in as a string (sapce delimited) concatenation of all sample.csv file
    # Ex.  sample1.csv sample2.csv sample3.csv 
    # which will then be split by the script to read each csv file
    combine_reports.py -o "report!{params.file_extension}.csv" -p "!{reports_file}" -i "!{params.reads}" -y report_mqc.yaml \
                       -q "!{primer_stats}" -z primer_stats_mqc.yaml -l !{params.primer} \
                       -r "!{read_length}" -x read_length_mqc.yaml

    '''

}