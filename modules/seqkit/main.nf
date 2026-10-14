process seqFilter {
        publishDir "${params.final_outdir}/seqFilter", mode: 'copy'
        input:
                path input_files
                val pattern
        output:

                path "${input_files.baseName}_filtered.fasta"
        script:
        """
ml seqkit/2.8.2

seqkit grep -nrp "${pattern}" ${input_files} > ${input_files.baseName}_filtered.fasta
        """
}
