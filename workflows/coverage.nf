#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MAKE_WHO_TARGETS {
    publishDir "${params.results_dir}/qc", mode: 'copy'

    input:
    path annotation_table
    path repetitive_regions
    path reference

    output:
    path "who_targets.bed", emit: targets

    script:
    def padding = params.promoter_padding ?: 200
    def mask_arg = repetitive_regions.name != 'NO_FILE' ? "--repetitive-regions ${repetitive_regions}" : ''
    """
    if [ ! -f ${reference}.fai ]; then
        samtools faidx ${reference}
    fi

    python3 ${baseDir}/scripts/make_who_targets.py \\
        --annotation-table ${annotation_table} \\
        --output who_targets.bed \\
        --promoter-padding ${padding} \\
        --reference-fai ${reference}.fai \\
        ${mask_arg}
    """

    stub:
    """
    touch who_targets.bed
    """
}

process MOSDEPTH {
    tag "${sample_id}"
    publishDir "${params.results_dir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path targets

    output:
    tuple val(sample_id), path("${sample_id}.coverage.json"), emit: coverage
    path "${sample_id}.regions.bed.gz",                        emit: regions,    optional: true
    path "${sample_id}.thresholds.bed.gz",                     emit: thresholds, optional: true

    script:
    def min_depth = params.coverage_min_depth ?: 10
    def breadth   = params.coverage_min_breadth ?: 0.95
    def threshold_list = ([min_depth as int, 30] as Set).sort().join(',')
    """
    if [ ! -f ${bam}.bai ] && [ ! -f ${bam.baseName}.bai ]; then
        samtools index ${bam}
    fi

    mosdepth \\
        --by ${targets} \\
        --thresholds ${threshold_list} \\
        --no-per-base \\
        --threads ${task.cpus} \\
        ${sample_id} \\
        ${bam}

    python3 ${baseDir}/scripts/summarize_coverage.py \\
        --sample-id ${sample_id} \\
        --targets ${targets} \\
        --regions ${sample_id}.regions.bed.gz \\
        --thresholds ${sample_id}.thresholds.bed.gz \\
        --min-depth ${min_depth} \\
        --min-breadth ${breadth} \\
        --output ${sample_id}.coverage.json
    """

    stub:
    """
    echo '{"sample_id": "${sample_id}", "coverage_assessed": false}' > ${sample_id}.coverage.json
    """
}

process COVERAGE_NOT_ASSESSED {
    tag "${sample_id}"
    publishDir "${params.results_dir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.coverage.json"), emit: coverage

    script:
    """
    python3 ${baseDir}/scripts/summarize_coverage.py \\
        --sample-id ${sample_id} \\
        --not-assessed "VCF-only input: no alignment available, coverage not assessed" \\
        --output ${sample_id}.coverage.json
    """

    stub:
    """
    echo '{"sample_id": "${sample_id}", "coverage_assessed": false}' > ${sample_id}.coverage.json
    """
}

workflow COVERAGE {
    take:
    bam_ch         
    vcf_only_ch    

    main:
    annotation_table_ch = Channel.fromPath(params.annotation_table, checkIfExists: true).first()
    repetitive_ch = Channel
        .fromPath(params.repetitive_regions, checkIfExists: false)
        .ifEmpty(file("${baseDir}/data/NO_FILE"))
        .first()
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true).first()

    targets = MAKE_WHO_TARGETS(annotation_table_ch, repetitive_ch, reference_ch)

    from_bam = MOSDEPTH(bam_ch, targets.targets)
    from_vcf = COVERAGE_NOT_ASSESSED(vcf_only_ch)

    all_coverage = from_bam.coverage.mix(from_vcf.coverage)

    emit:
    coverage = all_coverage
    targets  = targets.targets
}
