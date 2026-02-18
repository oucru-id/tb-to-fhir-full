#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CREATE_FHIR {
    publishDir "${params.results_dir}/fhir", mode: 'copy'

    input:
    path(annotation)
    path(lineage_files)

    output:
    path "*.fhir.json", emit: fhir_output
    path "versions.yml", emit: versions

    script:
    def sample_id = annotation.simpleName.replaceAll(/\.annotated_variants$/, '')
    """
    if [[ "${annotation}" == *.gz ]]; then
        gunzip -c ${annotation} > ${sample_id}.vcf
    else
        cp ${annotation} ${sample_id}.vcf
    fi

    mkdir -p lineage_data
    
    for file in ${lineage_files}; do
        if [ -f "\$file" ]; then
            cp "\$file" lineage_data/
        fi
    done

    python3 $baseDir/scripts/annotated_to_fhir.py \\
        --input ${sample_id}.vcf \\
        --output ${sample_id}.fhir.json \\
        --lineage_dir lineage_data/

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

    main:
    
    CREATE_FHIR(annotated_ch, lineage_ch)

    emit:
    fhir_output = CREATE_FHIR.out.fhir_output
    versions    = CREATE_FHIR.out.versions
}
