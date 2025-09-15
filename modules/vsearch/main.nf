process quality_filtering {
    publishDir "${params.final_outdir}/${sample}/temp", mode: 'copy', pattern: "*.fasta"
    tag "${sample}"
    // debug true

    input:
    tuple val(sample), path (fastq)

    output:
    tuple val(sample), path ("${sample}.fasta"), emit: fasta, optional:true
    path ("${sample}_qfilter_log.csv"), emit: log_csv, optional: true
    path "versions.yml"         , emit: versions, optional: true

    shell:
    """
    vsearch --fastx_filter !{fastq} --fastq_maxee 1 --fastaout !{sample}.fasta \
                            --log !{sample}_qfilter.log
    parse_qfilter_log.py -p !{sample}_qfilter.log -s !{sample} -o !{sample}_qfilter_log.csv

    #remove space between seq_id and =adapter 
    #choose not to use linux sed , because it could be very slow
    remove_space.py -f !{sample}.fasta

    echo "!{task.process}:" > versions.yml
    echo "  vsearch: \$(vsearch 2>&1 | head -n 1 | cut -d ' ' -f2 | cut -d '_' -f1 | sed 's/^v//')" >> versions.yml
    """

}

process dereplication {
    publishDir "${params.final_outdir}/${sample}/temp", mode: 'copy', pattern: "*.fasta"
    tag "${sample}"

    input:
    tuple val(sample), path (fasta)

    output:
    tuple val(sample), path ("${sample}.unique.fasta"), emit: fasta, optional:true
    path ("${sample}_dereplication_log.csv"), emit: log_csv, optional: true


    shell:
    '''
    vsearch --derep_fulllength !{fasta} --output !{sample}.unique.fasta --sizeout --relabel_keep \
            --log !{sample}_dereplication.log
    parse_derep_log.py -p !{sample}_dereplication.log -s !{sample} -o !{sample}_dereplication_log.csv

    '''

}

process denoising {
    publishDir "${params.final_outdir}/${sample}", mode: 'copy', pattern: "*.fasta"
    tag "${sample}"
    // debug true
    cpus = "${params.mincpus}"

    input:
    tuple val(sample), path (fasta), path(ch_primer_file)

    output:
    tuple val(sample), path ("${sample}.final.unique.fasta"), emit:unique, optional:true
    path ("${sample}_denoise_log.csv"), emit: log_csv, optional: true

    shell:
    '''
    run_cluster_unoise.sh !{fasta} !{ch_primer_file} !{sample}.unique.unoise.fasta \
                            !{params.denoising_minsize} !{params.denoising_alpha}
    
    #remove potential chimeras
    if [  -s "!{sample}.unique.unoise.fasta" ]; then
        vsearch --uchime3_denovo !{sample}.unique.unoise.fasta --nonchimeras !{sample}.final.unique.fasta
    fi

    #generate denoise_log.csv file
    # Check if the fasta file exists and is not empty
    if [ ! -s "!{fasta}" ]; then
        seq_count_before="n/a"
        seq_count_removed="n/a"
    else
        seq_count_before=$(grep -c '>' "!{fasta}")
    fi

    # Check if the final unique fasta file exists and is not empty
    if [ ! -s "!{sample}.final.unique.fasta" ]; then
        seq_count_after=0
    else
        seq_count_after=$(grep -c '>' "!{sample}.final.unique.fasta")
    fi

    # Calculate seq_count_removed if seq_count_before is a number
    if [[ "$seq_count_before" =~ ^[0-9]+$ ]]; then
        seq_count_removed=$((seq_count_before - seq_count_after))
    fi
    
    echo "Sample name, Total Reads, Removed Reads" > !{sample}_denoise_log.csv
    echo "!{sample},$seq_count_before,$seq_count_removed" >> !{sample}_denoise_log.csv

    '''

}

process search_exact {
    publishDir "${params.final_outdir}/${sample}/temp", mode: 'copy'
    tag "${sample}"
    cpus = "${params.medcpus}"

    input:
    tuple val(sample), path (output_fasta), path (final_unique_fasta)

    output:
    tuple val(sample), path ("${sample}.output.match.final.txt")

    shell:
    '''
    vsearch --search_exact !{output_fasta} -db !{final_unique_fasta} \
            --userfields target+query --userout  !{sample}.output.match.final.txt

    '''

}