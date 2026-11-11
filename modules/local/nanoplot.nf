process NANOPLOT {
    tag "${meta.id}"
    label 'process_low'
    conda "bioconda::nanoplot=1.48.0"
    container "quay.io/biocontainers/nanoplot:1.48.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_NanoStats.txt"),        emit: stats
    tuple val(meta), path("${meta.id}_NanoPlot-report.html"), emit: report
    path "versions.yml",                                      emit: versions

    script:
    """
    NanoPlot --fastq ${reads} --threads ${task.cpus} --outdir . --prefix '${meta.id}_' --no_static

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version 2>&1 | sed 's/^NanoPlot //')
    END_VERSIONS
    """

    stub:
    """
    touch '${meta.id}_NanoStats.txt' '${meta.id}_NanoPlot-report.html'
    printf '"%s":\\n    nanoplot: stub\\n' "${task.process}" > versions.yml
    """
}
