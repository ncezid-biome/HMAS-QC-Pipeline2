process multiqc {
    publishDir "${params.outdir}", mode: 'copy'
    tag "multiqc"

    input:
    path ('*')
    path(multiqc_config)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , optional:true, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    """
    multiqc \\
        --force \\
        $config \\
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