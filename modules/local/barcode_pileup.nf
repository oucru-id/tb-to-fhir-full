process BARCODE_PILEUP {
    tag "${meta.id}"
    label 'process_single'
    conda "bioconda::bcftools=1.19 bioconda::htslib=1.19.1"
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_1"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path barcode_bed

    output:
    tuple val(meta), path("${meta.id}.barcode_pileup.tsv"), emit: tsv
    path "versions.yml",                                    emit: versions

    script:
    def min_base_quality = meta.platform == 'ont' ? 7 : 13
    """
    bcftools mpileup -f ${fasta} -R ${barcode_bed} -a FORMAT/AD -q 20 -Q ${min_base_quality} -B -I -d 100000 -Ou ${bam} \\
        | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%AD]\\n' > ${meta.id}.barcode_pileup.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.barcode_pileup.tsv
    printf '"%s":\\n    bcftools: stub\\n' "${task.process}" > versions.yml
    """
}
