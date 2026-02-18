#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CONVERT_DEEPLEX {
    publishDir "${params.results_dir}/fhir_deeplex", mode: 'copy'
    
    input:
    path(deeplex_file)

    output:
    path "fhir_output/*.json", emit: fhir_files

    script:
    """
    mkdir -p fhir_output
    python3 $baseDir/scripts/xlsx_json_converter.py \\
        --input "${deeplex_file}" \\
        --output_dir fhir_output
    """
}

process MERGE_CLINICAL_DEEPLEX {
    publishDir "${params.results_dir}/fhir_deeplex_merged", mode: 'copy'

    input:
    path(fhir_bundle)
    path(clinical_metadata)

    output:
    path "*.merged.fhir.json", emit: merged_fhir

    script:
    def prefix = fhir_bundle.simpleName
    """
    python3 $baseDir/scripts/merge_clinical_deeplex.py \\
        --input "${fhir_bundle}" \\
        --output "${prefix}.merged.fhir.json" \\
        --clinical_metadata "${clinical_metadata}"
    """
}

workflow DEEPLEX {
    take:
    deeplex_ch
    clinical_metadata_ch

    main:
    fhir_raw = CONVERT_DEEPLEX(deeplex_ch)
    fhir_flattened = fhir_raw.fhir_files.flatten()
    merged_fhir = MERGE_CLINICAL_DEEPLEX(
        fhir_flattened,
        clinical_metadata_ch
    )

    emit:
    fhir_output = merged_fhir.merged_fhir
}