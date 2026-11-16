process CLASSIFY_LINEAGE {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12 conda-forge::openpyxl=3.1.5"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(vcf), path(pileup)
    path fasta
    path barcode_bed

    output:
    tuple val(meta), path("${meta.id}.lineage.json"), emit: json
    path "versions.yml",                              emit: versions

    script:
    def pileup_arg = pileup ? "--pileup ${pileup}" : ''
    """
    lineage_classifier.py \\
        --vcf ${vcf} \\
        ${pileup_arg} \\
        --barcode ${barcode_bed} \\
        --reference ${fasta} \\
        --min_depth ${params.filter_min_depth} \\
        --sample_id '${meta.id}' \\
        --output ${meta.id}.lineage.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        lineage_barcode: ${barcode_bed.name}
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id": "${meta.id}", "lineage": "unknown"}' > ${meta.id}.lineage.json
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
