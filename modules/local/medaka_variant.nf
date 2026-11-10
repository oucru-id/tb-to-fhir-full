process MEDAKA_VARIANT {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::medaka=2.1.0"
    container "quay.io/biocontainers/medaka:2.1.0--py38ha0c3a46_0"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta

    output:
    tuple val(meta), path("${meta.id}.variants.vcf.gz"), emit: vcf
    path "versions.yml",                                 emit: versions

    script:
    """
    cp -L ${fasta} reference.fasta
    samtools faidx reference.fasta
    export CUDA_VISIBLE_DEVICES=""

    medaka inference ${bam} consensus_probs.hdf \\
        --model ${params.medaka_model} \\
        --batch_size ${params.medaka_batch_size} \\
        --threads ${task.cpus}
    medaka vcf consensus_probs.hdf reference.fasta medaka.vcf
    bcftools sort medaka.vcf -o medaka.sorted.vcf
    medaka tools annotate medaka.sorted.vcf reference.fasta ${bam} medaka.annotated.vcf

    if [ ! -s medaka.annotated.vcf ]; then
        echo "ERROR: medaka produced no annotated VCF" >&2
        exit 1
    fi
    bgzip -c medaka.annotated.vcf > ${meta.id}.variants.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$(medaka --version 2>&1 | sed 's/medaka //')
        medaka_model: ${params.medaka_model}
    END_VERSIONS
    """

    stub:
    """
    echo | gzip > ${meta.id}.variants.vcf.gz
    printf '"%s":\\n    medaka: stub\\n' "${task.process}" > versions.yml
    """
}
