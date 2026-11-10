process FASTQC {
    tag "${meta.id}"
    label 'process_low'
    conda "bioconda::fastqc=0.12.1"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.zip"),  emit: zip
    tuple val(meta), path("*_fastqc.html"), emit: html
    path "versions.yml",                    emit: versions

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/FastQC v//')
    END_VERSIONS
    """

    stub:
    """
    for r in ${reads}; do
        b=\$(basename "\$r" .fastq.gz)
        touch "\${b}_fastqc.zip" "\${b}_fastqc.html"
    done
    printf '"%s":\\n    fastqc: stub\\n' "${task.process}" > versions.yml
    """
}
