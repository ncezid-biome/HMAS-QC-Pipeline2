process cutadapt {
    publishDir "${params.final_outdir}/${sample}", mode: 'copy', pattern: "cutadapt/*.json"
    publishDir "${params.final_outdir}/${sample}", mode: 'copy', pattern: "cutadapt/*.trimmed.*.fastq", enabled: params.save_trimmed
    publishDir "${params.final_outdir}/${sample}", mode: 'copy', pattern: "cutadapt/*.matched.*.fastq", enabled: params.save_trimmed
    publishDir "${params.final_outdir}/${sample}", mode: 'copy', pattern: "cutadapt/*.discarded.*.fastq", enabled: params.save_trimmed
    publishDir "${params.final_outdir}/${sample}", mode: 'copy', pattern: "cutadapt/*.adapter_filter.csv"
    tag "${sample}"
    cpus = "${params.maxcpus}"
    memory = "${params.medmems}"
    errorStrategy 'retry'
    maxRetries 3

    maxForks = "${params.maxcutadapts}"

    input:
    tuple val(sample), path(reads), path(ch_primer_file)

    output:
    tuple val(sample), path("cutadapt/${sample}.matched.1.fastq"), path("cutadapt/${sample}.matched.2.fastq"), emit: cutadapt_fastq, optional: true
    tuple val(sample), path("cutadapt/${sample}.cutadapt.json"), emit: cutadapt_json, optional: true
    path("cutadapt/${sample}.adapter_filter.csv"), emit: adapter_filter_csv, optional: true
    path("cutadapt/${sample}.trimmed.1.fastq"), emit: trimmed_r1, optional: true
    path("cutadapt/${sample}.trimmed.2.fastq"), emit: trimmed_r2, optional: true
    path("cutadapt/${sample}.discarded.1.fastq"), emit: discarded_r1, optional: true
    path("cutadapt/${sample}.discarded.2.fastq"), emit: discarded_r2, optional: true

    shell:
    '''
    mkdir -p cutadapt
    run_cutadapt.py -f !{reads[0]} -r !{reads[1]} \
                    -o cutadapt -s !{sample} -p !{ch_primer_file} \
                    -e !{params.cutadapt_maxerror} -l !{params.cutadapt_minlength} \
                    -t !{params.cutadapt_thread} -b !{params.cutadapt_long}

    filter_mismatched_adapters.py \
                    -r1 cutadapt/!{sample}.trimmed.1.fastq \
                    -r2 cutadapt/!{sample}.trimmed.2.fastq \
                    -s !{sample} \
                    -o cutadapt
    '''
}
