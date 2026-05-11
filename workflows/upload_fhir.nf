#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process UPLOAD_TO_FHIR {
    publishDir "${params.results_dir}/fhir_upload", mode: 'copy'
    debug true

    input:
    path(fhir_file)

    output:
    path "${fhir_file.baseName}.upload.json", emit: upload_result

    script:
    def token_file_arg = "--token_file '${projectDir}/data/access_token.json'"
    def api_key_arg    = params.api_key ? "--api_key '${params.api_key}'" : ""
    """
    python3 ${projectDir}/scripts/upload_fhir.py \\
        --fhir_file ${fhir_file} \\
        --fhir_server_url '${params.fhir_server_url}' \\
        ${token_file_arg} \\
        ${api_key_arg}
    """
}

workflow UPLOAD_FHIR {
    take:
    validated_fhir_files

    main:
    UPLOAD_TO_FHIR(validated_fhir_files)
    
    emit:
    results = UPLOAD_TO_FHIR.out.upload_result
}