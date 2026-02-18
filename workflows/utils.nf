nextflow.enable.dsl = 2

process VERSIONS {
    publishDir "${params.results_dir}", mode: 'copy'

    output:
    path "software_versions.yml"

    script:
    """
    echo "pipeline:" > software_versions.yml
    echo "  name: tb_mutation_analysis" >> software_versions.yml
    echo "  version: ${params.version}" >> software_versions.yml
    echo "  nextflow: $nextflow.version" >> software_versions.yml
    
    echo "databases:" >> software_versions.yml
    echo "  reference: \$(basename ${params.reference})" >> software_versions.yml
    echo "  mutation_db: \$(basename ${params.mutation_db})" >> software_versions.yml
    echo "  repetitive_regions: \$(basename ${params.repetitive_regions})" >> software_versions.yml
    
    echo "processing_settings:" >> software_versions.yml
    echo "  illumina:" >> software_versions.yml
    echo "    trimmomatic: 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36'" >> software_versions.yml
    echo "    gatk_min_base_qual: 10" >> software_versions.yml
    echo "    gatk_min_confidence: 10.0" >> software_versions.yml
    echo "    filter_min_depth: 5" >> software_versions.yml
    echo "    filter_min_gq: 20" >> software_versions.yml
    echo "  nanopore:" >> software_versions.yml
    echo "    chopper_min_q: ${params.nanopore_min_q}" >> software_versions.yml
    echo "    chopper_min_l: ${params.nanopore_min_l}" >> software_versions.yml
    echo "    medaka_model: r941_e81_sup_variant_g514" >> software_versions.yml
    echo "    filter_min_depth: 5" >> software_versions.yml
    echo "    filter_min_gq: 20" >> software_versions.yml
    
    export BASE_DIR="${baseDir}"
    python3 $baseDir/scripts/get_versions.py >> software_versions.yml
    """
}