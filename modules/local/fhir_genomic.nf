
process FHIR_GENOMIC {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12 conda-forge::openpyxl=3.1.5"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(annotated_vcf), path(lineage_json), path(coverage_json)
    path gene_map

    output:
    tuple val(meta), path("${meta.id}.genomic.fhir.json"), emit: json
    path "versions.yml",                                   emit: versions

    script:
    def args = []
    if (gene_map.name != 'NO_FILE')        args << "--gene_map ${gene_map}"
    if (params.performer_org_id)           args << "--performer_org_id '${params.performer_org_id}'"
    if (params.performer_org_name)         args << "--performer_org_name '${params.performer_org_name}'"
    if (params.identifier_system_base)     args << "--identifier_system_base '${params.identifier_system_base}'"
    if (params.legacy_annotation_table)    args << "--legacy-annotation-table"
    """
    tbfhir_genomic.py \\
        --input ${annotated_vcf} \\
        --sample_id '${meta.id}' \\
        --lineage_json ${lineage_json} \\
        --coverage_json ${coverage_json} \\
        --output ${meta.id}.genomic.fhir.json \\
        --fhir_profile ${params.fhir_profile} \\
        --pipeline_version '${workflow.manifest.version}' \\
        --filter_min_depth ${params.filter_min_depth} \\
        --filter_min_quality ${params.filter_min_quality} \\
        --coverage_min_depth ${params.coverage_min_depth} \\
        --coverage_min_breadth ${params.coverage_min_breadth} \\
        ${args.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        fhir_profile: ${params.fhir_profile}
    END_VERSIONS
    """

    stub:
    """
    echo '{"resourceType": "Bundle", "type": "collection", "entry": []}' > ${meta.id}.genomic.fhir.json
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
