process GATK4_HAPLOTYPECALLER {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::gatk4=4.4.0.0"
    container "quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("${meta.id}.variants.vcf.gz"), emit: vcf
    path "versions.yml",                                 emit: versions

    script:
    def mem = (task.memory.toMega() * 0.8).intValue()
    """
    gatk --java-options "-Xmx${mem}M" HaplotypeCaller \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${meta.id}.variants.vcf.gz \\
        --create-output-variant-index true \\
        --native-pair-hmm-threads ${task.cpus} \\
        --min-base-quality-score 10 \\
        --standard-min-confidence-threshold-for-calling 10.0 \\
        --annotation MappingQuality \\
        --annotation QualByDepth \\
        --annotation FisherStrand \\
        --output-mode EMIT_VARIANTS_ONLY \\
        --dont-use-soft-clipped-bases \\
        --assembly-region-padding 100 \\
        --max-reads-per-alignment-start 0 \\
        --min-pruning 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p')
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.variants.vcf.gz
    printf '"%s":\\n    gatk4: stub\\n' "${task.process}" > versions.yml
    """
}
