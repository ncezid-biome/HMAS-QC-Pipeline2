process combine_reports {
    publishDir "${params.final_outdir}", mode: 'copy', pattern: "*.csv"
    tag "combine reports"
    // debug true

    input:
    path (reports_file)
    path (primer_stats)
    path (read_length)
    path (ch_primer_file)
    val (sample_ids)

    output:
    path ("report*.csv"), emit: csv, optional:true
    path ("report_mqc.yaml"), emit: report_mqc, optional:true
    path ("primer_stats_mqc.yaml"), emit: primer_stats_mqc, optional:true
    path ("read_length_mqc.yaml"), emit: read_length_mqc, optional:true
    path ("versions.yml"), emit: versions, optional: true

    shell:
    '''
    # use inspect() to prepare a stringified representation of a list, and pass that to a python script
    # in python script, use ast.literal_eval() to safely evaluate a string containing a Python literal (a list in this case)
    combine_reports.py -o "report!{params.file_extension}.csv" -p "!{reports_file}" -i "!{sample_ids.inspect()}" -y report_mqc.yaml \
                       -q "!{primer_stats}" -z primer_stats_mqc.yaml -l !{ch_primer_file} \
                       -r "!{read_length}" -x read_length_mqc.yaml

    echo "!{task.process}:" > versions.yml
    echo "  python: $(python --version 2>&1 | awk '{print $2}')" >> versions.yml

    '''

}