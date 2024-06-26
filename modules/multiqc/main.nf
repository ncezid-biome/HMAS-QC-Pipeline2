process multiqc {
    publishDir "${params.outdir}", mode: 'copy'
    tag "multiqc"

    input:
    path ('*')
    path (config_yaml)
 
    output:
    path "multiqc_report.html", emit:multiqc_report, optional:true
    path "multiqc_data", emit: multiqc_data, optional:true
 
    script:
    """
    multiqc -c $config_yaml .
    """
}