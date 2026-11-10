process GATK4_CREATESEQUENCEDICTIONARY {
    tag "${fasta.name}"
    label 'process_single'
    conda "bioconda::gatk4=4.4.0.0"
    container "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"

    input:
    path fasta

    output:
    path "${fasta.baseName}.dict", emit: dict
    path "versions.yml",           emit: versions

    script:
    def mem = (task.memory.toMega() * 0.8).intValue()
    """
    gatk --java-options "-Xmx${mem}M" CreateSequenceDictionary -R ${fasta} -O ${fasta.baseName}.dict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta.baseName}.dict
    printf '"%s":\\n    gatk4: stub\\n' "${task.process}" > versions.yml
    """
}
