process ANNOTATE_VARIANTS {
    tag "${meta.id}"
    label 'process_single'
    conda "bioconda::bcftools=1.19 bioconda::htslib=1.19.1"
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_1"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path annotation_table
    path annotation_table_index
    path annotation_header

    output:
    tuple val(meta), path("${meta.id}.annotated_variants.vcf.gz"), emit: vcf
    path "versions.yml",                                           emit: versions

    script:
    """
    bcftools annotate \\
        -a ${annotation_table} \\
        -h ${annotation_header} \\
        -c CHROM,POS,REF,ALT,INFO/GENE,INFO/DRUG,INFO/EFFECT,INFO/WHO_CLASSIFICATION,INFO/VARIANT_ID,INFO/GENOME_POSITION \\
        -O z \\
        -o ${meta.id}.annotated_variants.vcf.gz \\
        ${vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
        who_annotation_table: ${annotation_table.name}
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.annotated_variants.vcf.gz
    printf '"%s":\\n    bcftools: stub\\n' "${task.process}" > versions.yml
    """
}
