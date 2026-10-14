process pyReport {
	publishDir "${params.final_outdir}/pyReport", mode: 'copy'
	input:
		path input_files
	output:
		path "${input_files.baseName}_report.tsv"
	script:
	"""
	ml python/3.7.4
	create_primer_report.py --input ${input_files} --output ${input_files.baseName}_report.tsv
	"""
}
