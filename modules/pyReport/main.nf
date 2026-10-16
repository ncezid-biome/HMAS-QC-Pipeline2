process pyReport {
	publishDir "${params.final_outdir}/pyReport", mode: 'copy'
	tag "${input_files}"
	input:
		path input_files
	output:
		path "${input_files.baseName}_report.tsv"
	script:
	"""
	create_primer_report.py --input ${input_files} --output ${input_files.baseName}_report.tsv
	"""
}
