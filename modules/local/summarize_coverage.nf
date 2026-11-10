process SUMMARIZE_COVERAGE {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12 conda-forge::openpyxl=3.1.5"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(regions), path(thresholds)
    path targets

    output:
    tuple val(meta), path("${meta.id}.coverage.json"), emit: json
    path "versions.yml",                               emit: versions

    script:
    """
    summarize_coverage.py \\
        --sample-id '${meta.id}' \\
        --targets ${targets} \\
        --regions ${regions} \\
        --thresholds ${thresholds} \\
        --min-depth ${params.coverage_min_depth} \\
        --min-breadth ${params.coverage_min_breadth} \\
        --output ${meta.id}.coverage.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id": "${meta.id}", "coverage_assessed": false, "reason": "stub run", "genes": {}, "drugs": {}}' > ${meta.id}.coverage.json
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}

process COVERAGE_NOT_ASSESSED {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12 conda-forge::openpyxl=3.1.5"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    val meta

    output:
    tuple val(meta), path("${meta.id}.coverage.json"), emit: json
    path "versions.yml",                               emit: versions

    script:
    """
    summarize_coverage.py \\
        --sample-id '${meta.id}' \\
        --not-assessed "VCF-only input: no alignment available, coverage not assessed" \\
        --min-depth ${params.coverage_min_depth} \\
        --min-breadth ${params.coverage_min_breadth} \\
        --output ${meta.id}.coverage.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id": "${meta.id}", "coverage_assessed": false, "reason": "stub run", "genes": {}, "drugs": {}}' > ${meta.id}.coverage.json
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
