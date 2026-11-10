
process FILTER_VARIANTS {
    tag "${meta.id}"
    label 'process_single'
    conda "bioconda::bcftools=1.19 bioconda::htslib=1.19.1"
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_1"

    input:
    tuple val(meta), path(vcf)
    path repetitive_regions

    output:
    tuple val(meta), path("${meta.id}.filtered_variants.vcf.gz"), path("${meta.id}.filtered_variants.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def dp = params.filter_min_depth
    def q  = params.filter_min_quality
    def out = "${meta.id}.filtered_variants.vcf.gz"
    def mask = repetitive_regions.name == 'NO_FILE'
        ? "echo 'WARNING: no repetitive-regions mask (params.repetitive_regions = NO_FILE)' >&2\n    cp ${vcf} step1.vcf.gz"
        : "bcftools view -T ^${repetitive_regions} ${vcf} -O z -o step1.vcf.gz"
    def filters
    if (meta.platform == 'illumina') {
        filters = """
    bcftools view -v snps,indels -i 'FORMAT/DP>=0 && (GT="0/0" || GT="1/1" || GT="0/1")' -O z -o step2.vcf.gz step1.vcf.gz
    bcftools view -i "FORMAT/DP>=${dp}" -O z -o step3.vcf.gz step2.vcf.gz
    if zcat step3.vcf.gz | grep -q "^##FORMAT=<ID=GQ"; then
        bcftools view -i "FORMAT/GQ>=${q}" -O v step3.vcf.gz | bcftools sort -O z -o ${out}
    else
        bcftools sort -O z -o ${out} step3.vcf.gz
    fi"""
    } else if (meta.platform == 'ont') {
        filters = """
    bcftools view -v snps,indels -i 'DP>=0 && (GT="0" || GT="1" || GT="0/1")' -O z -o step2.vcf.gz step1.vcf.gz
    bcftools view -i "DP>=${dp}" -O z -o step3.vcf.gz step2.vcf.gz
    bcftools view -i "GQ>=${q}" -O v step3.vcf.gz | bcftools sort -O z -o ${out}"""
    } else {
        // v1.6 used "FORMAT/GQ || INFO/GQ || QUAL" whenever either GQ header was
        // present; bcftools rejects an expression naming a tag the header lacks,
        // so the expression is now built from the tags that exist.
        filters = """
    bcftools view -v snps,indels -i 'FORMAT/DP>=0 && (GT="0/0" || GT="1/1" || GT="0/1" || GT="0" || GT="1")' -O z -o step2.vcf.gz step1.vcf.gz
    bcftools view -i "FORMAT/DP>=${dp}" -O z -o step3.vcf.gz step2.vcf.gz
    expr="QUAL>=${q}"
    if zcat step3.vcf.gz | grep -q "^##INFO=<ID=GQ,"; then expr="INFO/GQ>=${q} || \${expr}"; fi
    if zcat step3.vcf.gz | grep -q "^##FORMAT=<ID=GQ,"; then expr="FORMAT/GQ>=${q} || \${expr}"; fi
    bcftools view -i "\${expr}" -O v step3.vcf.gz | bcftools sort -O z -o ${out}"""
    }
    """
    ${mask}
    ${filters}
    tabix -p vcf ${out}
    rm -f step1.vcf.gz step2.vcf.gz step3.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.filtered_variants.vcf.gz
    touch ${meta.id}.filtered_variants.vcf.gz.tbi
    printf '"%s":\\n    bcftools: stub\\n' "${task.process}" > versions.yml
    """
}
