process blastQuery {
        publishDir "${params.final_outdir}/blastQuery", mode: 'copy'
	container 'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759'
	tag "${input_files}"
	 input:
                path input_files
                path db_dir
                val db_name

        output:

                path "${input_files.baseName}_query.tsv"
        script:
        """
blastn -query ${input_files} -db ${db_dir}/${db_name} -out ${input_files.baseName}_query.tsv -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

        """
}
