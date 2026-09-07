#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CREATE_FHIR {
    publishDir "${params.results_dir}/fhir", mode: 'copy'

    input:
    path(annotation)
    path(lineage_files)
    path(coverage_files)
    path(org_metadata)

    output:
    path "*.fhir.json", emit: fhir_output
    path "versions.yml", emit: versions

    script:
    def sample_id = annotation.simpleName.replaceAll(/\.annotated_variants$/, '')
    def min_depth = params.coverage_min_depth ?: 10
    def min_breadth = params.coverage_min_breadth ?: 0.95
    """
    if [[ "${annotation}" == *.gz ]]; then
        gunzip -c ${annotation} > ${sample_id}.vcf
    else
        cp ${annotation} ${sample_id}.vcf
    fi

    mkdir -p lineage_data coverage_data

    for file in *.lineage.json; do
        if [ -f "\$file" ]; then
            cp "\$file" lineage_data/
        fi
    done

    for file in *.coverage.json; do
        if [ -f "\$file" ]; then
            cp "\$file" coverage_data/
        fi
    done

    python3 $baseDir/scripts/annotated_to_fhir.py \\
        --input ${sample_id}.vcf \\
        --output ${sample_id}.fhir.json \\
        --lineage_dir lineage_data/ \\
        --coverage_dir coverage_data/ \\
        --pipeline_version '${params.version}' \\
        --coverage_min_depth ${min_depth} \\
        --coverage_min_breadth ${min_breadth} \\
        --organization_metadata ${org_metadata}

    cat <<-END_VERSIONS > versions.yml
    "fhir_converter":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def sample_id = annotation.simpleName.replaceAll(/\.annotated_variants$/, '')
    """
    touch ${sample_id}.fhir.json
    touch versions.yml
    """
}

workflow FHIR {
    take:
    annotated_ch
    lineage_ch
    coverage_ch
    org_metadata_ch

    main:

    CREATE_FHIR(annotated_ch, lineage_ch, coverage_ch, org_metadata_ch)

    emit:
    fhir_output = CREATE_FHIR.out.fhir_output
    versions    = CREATE_FHIR.out.versions
}
