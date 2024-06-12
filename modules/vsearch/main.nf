process quality_filtering {
    publishDir "${params.outdir}/${sample}/temp", mode: 'copy'
    tag "${sample}"

    input:
    tuple val(sample), path (fastq)

    output:
    tuple val(sample), path ("${sample}.fasta"), optional:true

    shell:
    '''
    vsearch --fastx_filter !{fastq} --fastq_maxee 1 --fastaout !{sample}.fasta

    #remove space between seq_id and =adapter 
    #choose not to use linux sed , because it could be very slow
    remove_space.py -f !{sample}.fasta
    '''

}

process dereplication {
    publishDir "${params.outdir}/${sample}/temp", mode: 'copy'
    tag "${sample}"

    input:
    tuple val(sample), path (fasta)

    output:
    tuple val(sample), path ("${sample}.unique.fasta")

    shell:
    '''
    vsearch --derep_fulllength !{fasta} --output !{sample}.unique.fasta --sizeout --relabel_keep

    '''

}

process denoising {
    publishDir "${params.outdir}/${sample}", mode: 'copy'
    tag "${sample}"
    debug true
    cpus = "${params.mincpus}"

    input:
    tuple val(sample), path (fasta)

    output:
    tuple val(sample), path ("${sample}.final.unique.fasta"), emit:unique, optional:true

    shell:
    '''
    vsearch --cluster_unoise !{fasta} --minsize !{params.denoising_minsize} \
                                      --unoise_alpha !{params.denoising_alpha} \
                                      --centroids !{sample}.unique.unoise.fasta
    
    #remove potential chimeras
    vsearch --uchime3_denovo !{sample}.unique.unoise.fasta --nonchimeras !{sample}.final.unique.fasta

    '''

}

process search_exact {
    publishDir "${params.outdir}/${sample}/temp", mode: 'copy'
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