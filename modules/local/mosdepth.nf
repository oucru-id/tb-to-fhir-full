process MOSDEPTH {
    tag "${meta.id}"
    label 'process_low'
    conda "bioconda::mosdepth=0.3.6"
    container "quay.io/biocontainers/mosdepth:0.3.6--hd299d5a_0"

    input:
    tuple val(meta), path(bam), path(bai)
    path targets

    output:
    tuple val(meta), path("${meta.id}.regions.bed.gz"), path("${meta.id}.thresholds.bed.gz"), emit: beds
    path "versions.yml", emit: versions

    script:
    def thresholds = ([params.coverage_min_depth as int, 30] as Set).sort().join(',')
    """
    mosdepth \\
        --by ${targets} \\
        --thresholds ${thresholds} \\
        --no-per-base \\
        --threads ${task.cpus} \\
        ${meta.id} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.regions.bed.gz
    echo | gzip > ${meta.id}.thresholds.bed.gz
    printf '"%s":\\n    mosdepth: stub\\n' "${task.process}" > versions.yml
    """
}
