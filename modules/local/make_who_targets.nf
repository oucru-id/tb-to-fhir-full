process MAKE_WHO_TARGETS {
    label 'process_single'
    conda "conda-forge::python=3.12 conda-forge::openpyxl=3.1.5"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    path annotation_table
    path repetitive_regions
    path fai

    output:
    path "who_targets.bed",   emit: bed
    path "who_positions.bed", emit: positions
    path "versions.yml",      emit: versions

    script:
    def mask = repetitive_regions.name != 'NO_FILE' ? "--repetitive-regions ${repetitive_regions}" : ''
    """
    make_who_targets.py \\
        --annotation-table ${annotation_table} \\
        --output who_targets.bed \\
        --positions-output who_positions.bed \\
        --promoter-padding ${params.promoter_padding} \\
        --reference-fai ${fai} \\
        ${mask}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'NC_000962.3\\t0\\t100\\trpoB\\tRifampicin\\n' > who_targets.bed
    printf 'NC_000962.3\\t50\\t51\\trpoB\\tRifampicin\\n' > who_positions.bed
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
