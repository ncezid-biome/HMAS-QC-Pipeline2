process seqFilter {
        publishDir "${params.final_outdir}/seqFilter", mode: 'copy'
        container 'community.wave.seqera.io/library/seqkit:2.13.0--05c0a96bf9fb2751'
	tag "${sample}"
	input:
                tuple val(sample), path(input_files)
                val pattern
        output:

                tuple val(sample), path("${sample}_filtered.fasta")
        script:
        """
seqkit grep -nrp "${pattern}" ${input_files} > ${sample}_filtered.fasta
        """
}
