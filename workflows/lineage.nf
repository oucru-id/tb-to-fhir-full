#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process extract_barcode_snps {
    
    input:
    tuple val(sample_id), path(vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}.barcode_variants.vcf")
    
    script:
    def barcode_bed = params.lineage_barcode
    """
    if [[ ${vcf} == *.gz ]]; then
        if [ ! -f ${vcf}.tbi ]; then
            tabix -p vcf ${vcf}
        fi
        
        bcftools view -R ${barcode_bed} ${vcf} > ${sample_id}.barcode_variants.vcf
    else
        bcftools view -R ${barcode_bed} ${vcf} > ${sample_id}.barcode_variants.vcf
    fi
    
    """
}

process classify_lineage {
    publishDir "${params.results_dir}/lineage", mode: 'copy'
    
    input:
    tuple val(sample_id), path(barcode_vcf)
    
    output:
    tuple val(sample_id), path("${sample_id}.lineage.json")
    
    script:
    """
    python3 $baseDir/scripts/lineage_classifier.py \
        --vcf ${barcode_vcf} \
        --barcode ${params.lineage_barcode} \
        --sample_id ${sample_id} \
        --output ${sample_id}.lineage.json
    """
}

workflow LINEAGE {
    take:
    filtered_vcfs
    
    main:
    barcode_snps = extract_barcode_snps(filtered_vcfs)
    lineage_results = classify_lineage(barcode_snps)
    
    emit:
    lineage_results = lineage_results  
}