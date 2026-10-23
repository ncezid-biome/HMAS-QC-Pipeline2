process blastQuery {
        publishDir "${params.final_outdir}/${sample}/blastQuery", mode: 'copy'
	container 'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759'
	tag "${sample}"
	 input:
                tuple val(sample), path(input_files)
                path db_dir
                val db_name

        output:

               tuple val(sample), path ("${sample}_query.tsv")
        script:
        """
blastn -query ${input_files} -db ${db_dir}/${db_name} -out ${sample}_query.tsv -evalue 1e-5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'

        """
}
