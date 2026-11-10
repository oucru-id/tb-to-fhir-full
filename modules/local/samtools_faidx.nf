process SAMTOOLS_FAIDX {
    tag "${fasta.name}"
    label 'process_single'
    conda "bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"

    input:
    path fasta

    output:
    path "${fasta}.fai", emit: fai
    path "versions.yml", emit: versions

    script:
    """
    samtools faidx ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.fai
    printf '"%s":\\n    samtools: stub\\n' "${task.process}" > versions.yml
    """
}
