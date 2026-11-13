process CHOPPER {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::chopper=0.12.0"
    container "quay.io/biocontainers/chopper:0.12.0--hcdda2d0_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed.fastq.gz"), emit: reads
    path "versions.yml",                                  emit: versions

    script:
    """
    zcat ${reads} \\
        | chopper -q ${params.nanopore_min_q} -l ${params.nanopore_min_l} --threads ${task.cpus} \\
        | gzip > ${meta.id}_trimmed.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | sed 's/chopper //')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}_trimmed.fastq.gz
    printf '"%s":\\n    chopper: stub\\n' "${task.process}" > versions.yml
    """
}
