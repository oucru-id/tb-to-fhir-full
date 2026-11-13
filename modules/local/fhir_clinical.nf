
process FHIR_CLINICAL {
    label 'process_single'
    conda "conda-forge::python=3.12 conda-forge::openpyxl=3.1.5"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    path patient_metadata
    path organization_metadata
    path practitioner_metadata
    path sample_ids

    output:
    path "clinical",                       emit: dir
    path "clinical/clinical_summary.tsv",  emit: summary
    path "versions.yml",                   emit: versions

    script:
    def args = []
    if (organization_metadata.name != 'NO_FILE') args << "--organization_metadata ${organization_metadata}"
    if (practitioner_metadata.name != 'NO_FILE') args << "--practitioner_metadata ${practitioner_metadata}"
    if (params.identifier_system_base)           args << "--identifier_system_base '${params.identifier_system_base}'"
    """
    clinical_to_fhir.py \\
        --patient_metadata ${patient_metadata} \\
        --sample_ids ${sample_ids} \\
        --fhir_profile ${params.fhir_profile} \\
        --outdir clinical \\
        ${args.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir clinical
    echo '{"resourceType": "Bundle", "type": "collection", "entry": []}' > clinical/shared.clinical.json
    echo '{}' > clinical/external_refs.json
    printf 'sample_id\\tfhir_id\\tstatus\\tflags\\n' > clinical/clinical_summary.tsv
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
