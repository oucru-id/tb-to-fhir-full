#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MULTIQC {
    publishDir "${params.results_dir}/qc", mode: 'copy'
    debug true

    input:
    path(qc_files)

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p fastqc_results
    cp -L ${qc_files} fastqc_results/

    multiqc fastqc_results \
        --force \
        --verbose \
        --outdir . \
        --no-ansi \
        --interactive \
        --zip-data-dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //')
    END_VERSIONS
    """
}

process CREATE_SAMPLE_REPORTS {
    publishDir "${params.results_dir}/reports", mode: 'copy'
    debug true

    input:
    path(annotated_vcf)
    path(lineage_files)

    output:
    path "*.summary_report.txt", emit: summary_reports
    path "versions.yml", emit: versions

    script:
    def sample_id = annotated_vcf.simpleName.replaceAll(/\.annotated_variants$/, '')
    """
    
    mkdir -p lineage_data
    
    for file in ${lineage_files}; do
        if [ -f "\$file" ]; then
            echo "Copying lineage file: \$file"
            cp "\$file" lineage_data/
        fi
    done

    python3 $baseDir/scripts/generate_sample_report.py \\
        --annotated_vcf ${annotated_vcf} \\
        --lineage_dir lineage_data/ \\
        --output_dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}

workflow GENERATE_REPORT {
    take:
    qc_files

    main:
    MULTIQC(qc_files)

    emit:
    report = MULTIQC.out.report
    versions = MULTIQC.out.versions
}

workflow GENERATE_SAMPLE_REPORTS {
    take:
    annotated_ch
    lineage_ch

    main:
    CREATE_SAMPLE_REPORTS(annotated_ch, lineage_ch)

    emit:
    summary_reports = CREATE_SAMPLE_REPORTS.out.summary_reports
    versions = CREATE_SAMPLE_REPORTS.out.versions
}