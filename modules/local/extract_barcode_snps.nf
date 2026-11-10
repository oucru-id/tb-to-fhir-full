process EXTRACT_BARCODE_SNPS {
    tag "${meta.id}"
    label 'process_single'
    conda "bioconda::bcftools=1.19 bioconda::htslib=1.19.1"
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_1"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path barcode_bed

    output:
    tuple val(meta), path("${meta.id}.barcode_variants.vcf"), emit: vcf
    path "versions.yml",                                      emit: versions

    script:
    """
    bcftools view -R ${barcode_bed} ${vcf} > ${meta.id}.barcode_variants.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.barcode_variants.vcf
    printf '"%s":\\n    bcftools: stub\\n' "${task.process}" > versions.yml
    """
}
