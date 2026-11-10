process SAMPLE_REPORT {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12 conda-forge::openpyxl=3.1.5"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(annotated_vcf), path(lineage_json)

    output:
    path "${meta.id}.summary_report.txt", emit: report
    path "versions.yml",                  emit: versions

    script:
    """
    generate_sample_report.py \\
        --annotated_vcf ${annotated_vcf} \\
        --sample_id '${meta.id}' \\
        --lineage_json ${lineage_json} \\
        --output_dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.summary_report.txt
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
