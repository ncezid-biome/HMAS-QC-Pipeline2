process seqFilter {
        publishDir "${params.final_outdir}/seqFilter", mode: 'copy'
        container 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'
	tag "${input_files}"
	input:
                path input_files
                val pattern
        output:

                path "${input_files.baseName}_filtered.fasta"
        script:
        """
seqkit grep -nrp "${pattern}" ${input_files} > ${input_files.baseName}_filtered.fasta
        """
}
