process combine_logs {
    // publishDir "${params.final_outdir}", mode: 'copy'
    tag "combine logs ${module_name}"
    // debug true

    input:
    path (logs_file)
    val (module_name)

    output:
    path ("${module_name}_mqc.out"), emit: log, optional: true

    shell:
    '''
    #!{logs_file} is passed in as a string (sapce delimited) concatenation of all sample log file
    # Ex.  sample1.log sample2.log sample3.log 
    # which will then be split by the script to read each csv file
    combine_logs.py -o "!{module_name}_mqc.out" -p "!{logs_file}" 

    '''

}