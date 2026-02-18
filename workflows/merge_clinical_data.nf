#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MERGE_CLINICAL_FHIR {
    publishDir "${params.results_dir}/fhir_merged", mode: 'copy'

    input:
    path(fhir_bundle)
    each path(clinical_metadata) 

    output:
    path "*.merged.fhir.json", emit: merged_fhir
    path "versions.yml", emit: versions

    script:
    def prefix = fhir_bundle.simpleName
    """
    
    python3 $baseDir/scripts/merge_clinical_fhir.py \\
        --input ${fhir_bundle} \\
        --output ${prefix}.merged.fhir.json \\
        --clinical_metadata ${clinical_metadata}

    cat <<-END_VERSIONS > versions.yml
    "merge_clinical_fhir":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}

workflow MERGE_CLINICAL_DATA {
    take:
    fhir_ch
    clinical_ch

    main:
    MERGE_CLINICAL_FHIR(fhir_ch, clinical_ch)

    emit:
    merged_fhir = MERGE_CLINICAL_FHIR.out.merged_fhir
    versions    = MERGE_CLINICAL_FHIR.out.versions
}