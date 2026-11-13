
process FHIR_UPLOAD {
    tag "${bundle.name}"
    label 'process_single'
    conda "conda-forge::python=3.12 conda-forge::openpyxl=3.1.5"
    container "quay.io/biocontainers/python:3.12.12"
    maxForks 1
    secret 'TBFHIR_CLIENT_ID'
    secret 'TBFHIR_CLIENT_SECRET'

    input:
    path bundle

    output:
    path "*.upload.json", emit: result
    path "versions.yml",  emit: versions

    script:
    template 'upload.sh'

    stub:
    """
    echo '{"status": "stub"}' > ${bundle.name.replace('.fhir.json', '')}.upload.json
    printf '"%s":\\n    upload: stub\\n' "${task.process}" > versions.yml
    """
}

process FHIR_UPLOAD_WITH_API_KEY {
    tag "${bundle.name}"
    label 'process_single'
    conda "conda-forge::python=3.12 conda-forge::openpyxl=3.1.5"
    container "quay.io/biocontainers/python:3.12.12"
    maxForks 1
    secret 'TBFHIR_CLIENT_ID'
    secret 'TBFHIR_CLIENT_SECRET'
    secret 'TBFHIR_API_KEY'

    input:
    path bundle

    output:
    path "*.upload.json", emit: result
    path "versions.yml",  emit: versions

    script:
    template 'upload.sh'

    stub:
    """
    echo '{"status": "stub"}' > ${bundle.name.replace('.fhir.json', '')}.upload.json
    printf '"%s":\\n    upload: stub\\n' "${task.process}" > versions.yml
    """
}
