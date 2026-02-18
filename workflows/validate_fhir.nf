#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process VALIDATE_FHIR {
    publishDir "${params.results_dir}/fhir_validated", mode: 'copy'
    debug true
    
    input:
    path(fhir_file)
    path(validator)
    
    output:
    path "${fhir_file.baseName}.validation.txt", emit: validation_report
    path "${fhir_file}", emit: validated_fhir
    
    script:
    """
    java -jar ${validator} \\
        ${fhir_file} \\
        -version 3.0.0 \\
        -ig hl7.fhir.uv.genomics-reporting#current \\
        > ${fhir_file.baseName}.validation.txt 2>&1 || echo "Validation warnings/errors: ${fhir_file.baseName}.validation.txt" >&2
    """
}

workflow VALIDATE {
    take:
    fhir_json_files
    
    main:
    validator_path = file(params.fhir_validator)
    VALIDATE_FHIR(fhir_json_files.flatten(), validator_path)
    
    emit:
    validation_report = VALIDATE_FHIR.out.validation_report
    validated_fhir = VALIDATE_FHIR.out.validated_fhir
}