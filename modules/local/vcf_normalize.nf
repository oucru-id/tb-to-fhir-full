process VCF_NORMALIZE {
    tag "${meta.id}"
    label 'process_single'
    conda "bioconda::bcftools=1.19 bioconda::htslib=1.19.1"
    container "quay.io/biocontainers/bcftools:1.19--h8b25389_1"

    input:
    tuple val(meta), path(vcf)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.normalized.vcf.gz"), emit: vcf
    path "versions.yml",                                   emit: versions

    script:
    """
    if [[ ${vcf} == *.gz ]]; then
        input_vcf=${vcf}
    else
        bgzip -c ${vcf} > input.vcf.gz
        input_vcf=input.vcf.gz
    fi

    ref_chr=\$(head -1 ${fai} | cut -f1)
    vcf_chr=\$(zcat \${input_vcf} | grep -v "^#" | head -1 | cut -f1)
    if [ -n "\${vcf_chr}" ] && [ "\${vcf_chr}" != "\${ref_chr}" ]; then
        printf '%s\\t%s\\n' "\${vcf_chr}" "\${ref_chr}" > chr_rename.txt
        bcftools annotate --rename-chrs chr_rename.txt \${input_vcf} -O z -o renamed.vcf.gz
        input_vcf=renamed.vcf.gz
    fi

    bcftools sort -O z -o sorted.vcf.gz \${input_vcf}
    tabix -p vcf sorted.vcf.gz
    bcftools norm -m -both -f ${fasta} --check-ref s -O z -o ${meta.id}.normalized.vcf.gz sorted.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.normalized.vcf.gz
    printf '"%s":\\n    bcftools: stub\\n' "${task.process}" > versions.yml
    """
}
