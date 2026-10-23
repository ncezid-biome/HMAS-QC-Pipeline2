process pyReport {
	publishDir "${params.final_outdir}/${sample}/pyReport", mode: 'copy'
	tag "${sample}"
	input:
		tuple val(sample), path (input_files)
	output:
		path "${sample}_report.tsv"
	script:
	"""
	create_primer_report.py --input ${input_files} --output ${sample}_report.tsv
	"""
}
