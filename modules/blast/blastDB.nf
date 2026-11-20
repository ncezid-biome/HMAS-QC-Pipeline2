process blastDB {
    publishDir "${params.final_outdir}", mode: 'copy'
    container 'community.wave.seqera.io/library/blast:2.17.0--d4fb881691596759'
    tag "blastDB"
    input:
        path hashed_fasta

    output:
        // Emit the actual database files
        path ("blast_db"), emit: db_files
        // Emit the path name for the blast database
        // val "${hashed_fasta}", emit: db_name

    script:
    """
    mkdir -p blast_db
    makeblastdb \\
        -in ${hashed_fasta} \\
        -dbtype nucl \\
        -out blast_db/reference_db \\
        -parse_seqids
    """
}
