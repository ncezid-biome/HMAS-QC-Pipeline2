process combine_reports {
    publishDir "${params.outdir}", mode: 'copy'
    tag "combine reports"
    // debug true

    input:
    path (reports_file)

    output:
    path ("report.csv"), emit: csv, optional:true
    path ("report_mqc.yaml"), emit: report_mqc, optional:true

    shell:
    '''
    #!{reports_file} is passed in as a string (sapce delimited) concatenation of all sample.csv file
    # Ex.  sample1.csv sample2.csv sample3.csv 
    # which will then be split by the script to read each csv file
    combine_reports.py -o report.csv -p "!{reports_file}" -i "!{params.reads}" -y report_mqc.yaml

    '''

}