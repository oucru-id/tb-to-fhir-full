process MINIMAP2_ALIGN {
    tag "${meta.id}"
    label 'process_high'
    conda "bioconda::minimap2=2.26 bioconda::samtools=1.17"
    container "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:7e6194c85b2f194e301c71cdda1c002a754e8cc1-0"

    input:
    tuple val(meta), path(reads)
    path fasta

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    path "versions.yml",                                                 emit: versions

    script:
    """
    minimap2 -x map-ont --secondary=no -L --MD \\
        -A 2 -B 4 -O 4,24 -E 2,1 \\
        -t ${task.cpus} -a ${fasta} ${reads} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai
    printf '"%s":\\n    minimap2: stub\\n' "${task.process}" > versions.yml
    """
}
