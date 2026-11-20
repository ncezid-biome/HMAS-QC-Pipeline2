process hash_reference {
    tag "hash_reference"
    publishDir "${params.final_outdir}/reference", mode: 'copy'
    input:
    path reference_fasta

    output:
    path "reference.hashed.fasta", emit: hashed_fasta
    path "id_to_md5.tsv",          emit: translation_table

    script:
    """
    hash_reference.py \\
        --input ${reference_fasta} \\
        --out-fasta reference.hashed.fasta \\
        --out-table id_to_md5.tsv
    """
}
