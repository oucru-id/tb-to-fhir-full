process BWAMEM2_INDEX {
    tag "${fasta.name}"
    label 'process_medium'
    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19.2"
    container "quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0"

    input:
    path fasta

    output:
    path "bwamem2",      emit: index
    path "versions.yml", emit: versions

    script:
    """
    mkdir bwamem2
    bwa-mem2 index -p bwamem2/${fasta} ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>/dev/null | tail -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir bwamem2
    touch bwamem2/${fasta}.0123 bwamem2/${fasta}.amb bwamem2/${fasta}.ann bwamem2/${fasta}.bwt.2bit.64 bwamem2/${fasta}.pac
    printf '"%s":\\n    bwa-mem2: stub\\n' "${task.process}" > versions.yml
    """
}
