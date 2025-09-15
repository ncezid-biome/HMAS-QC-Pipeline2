process multiqc {
    publishDir "${params.final_outdir}", mode: 'copy'
    tag "multiqc"
    container 'https://depot.galaxyproject.org/singularity/multiqc%3A1.31--pyhdfd78af_0'
    maxRetries 2


    input:
    path ('*')
    path(multiqc_config)
    path(all_versions)

    output:
    path "*.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , optional:true, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mqc_file_name=\$(basename "$params.final_outdir")
    cp $multiqc_config ${multiqc_config}.bak
    update_multiqc_config.py ${multiqc_config}.bak '${params.multiqc_header}' $all_versions '${params.primer}'

    if [ -n "$multiqc_config" ]; then
        config="--config ${multiqc_config}.bak"
    else
        config=""
    fi

    multiqc \\
        -n "\${mqc_file_name}.html" \\
        --force \\
        \$config \\
        .

    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    #END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_plots
    touch multiqc_report.html

    #cat <<-END_VERSIONS > versions.yml
    #"${task.process}":
    #    multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    #END_VERSIONS
    """

}
