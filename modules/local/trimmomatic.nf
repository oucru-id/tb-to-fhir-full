process TRIMMOMATIC {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::trimmomatic=0.39"
    container "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_{1,2}.trimmed.fastq.gz"), emit: reads
    path "versions.yml",                                        emit: versions

    script:

    def cmd = params.trimmomatic_jar ? "java -jar ${params.trimmomatic_jar}" : "trimmomatic"
    """
    ${cmd} PE -threads ${task.cpus} -phred33 \\
        ${reads[0]} ${reads[1]} \\
        ${meta.id}_1.trimmed.fastq.gz ${meta.id}_1.unpaired.fastq.gz \\
        ${meta.id}_2.trimmed.fastq.gz ${meta.id}_2.unpaired.fastq.gz \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(${cmd} -version 2>/dev/null | tail -1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_1.trimmed.fastq.gz ${meta.id}_2.trimmed.fastq.gz
    printf '"%s":\\n    trimmomatic: stub\\n' "${task.process}" > versions.yml
    """
}
