process MULTIQC {
    label 'process_low'
    conda "bioconda::multiqc=1.30"
    container "quay.io/biocontainers/multiqc:1.30--pyhdfd78af_0"

    input:
    path qc_files, stageAs: 'fastqc/*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data*",       emit: data
    path "versions.yml",        emit: versions

    script:
    """
    multiqc fastqc --force --outdir . --no-ansi --interactive --zip-data-dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_data
    printf '"%s":\\n    multiqc: stub\\n' "${task.process}" > versions.yml
    """
}
