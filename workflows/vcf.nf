#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process VCF_NORMALIZE {

    input:
    tuple val(sample_id), path(vcf_file)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}_normalized.vcf.gz")

    script:
    """
    if [[ ${vcf_file} == *.gz ]]; then
        input_vcf=${vcf_file}
    else
        bgzip -c ${vcf_file} > ${vcf_file}.gz
        input_vcf=${vcf_file}.gz
    fi

    ref_chr=\$(grep "^>" ${reference} | head -1 | sed 's/^>//' | cut -d' ' -f1)
    vcf_chr=\$(zcat \${input_vcf} | grep -v "^#" | head -1 | cut -f1)
    
    if [ "\${vcf_chr}" != "\${ref_chr}" ]; then
        echo "\${vcf_chr}\t\${ref_chr}" > chr_rename.txt
        bcftools annotate --rename-chrs chr_rename.txt \${input_vcf} -O z -o renamed.vcf.gz
        input_vcf=renamed.vcf.gz
    fi

    bcftools sort -O z -o sorted.vcf.gz \${input_vcf}
    tabix -p vcf sorted.vcf.gz

    bcftools norm \\
        -m -both \\
        -f ${reference} \\
        --check-ref s \\
        -O z \\
        -o ${sample_id}_normalized.vcf.gz \\
        sorted.vcf.gz

    tabix -p vcf ${sample_id}_normalized.vcf.gz
    """
}

process filter_variants {

    input:
    tuple val(sample_id), path(vcf_file)
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered_variants.vcf.gz")

    script:
    def repetitive_regions = "${params.repetitive_regions}"
    """
    # Exclude repetitive regions (pe/ppe genes, insertion sequences)
    if [ -f ${repetitive_regions} ]; then
        bcftools view -T ^${repetitive_regions} ${vcf_file} -O z -o step1.vcf.gz
    else
        echo "Warning: Repetitive regions file not found, skipping region filtering"
        cp ${vcf_file} step1.vcf.gz
    fi
    
    # Filter variants
    bcftools view \\
        -v snps,indels \\
        -i 'FORMAT/DP>=0 && (GT="0/0" || GT="1/1" || GT="0/1" || GT="0" || GT="1")' \\
        -O z \\
        -o step2.vcf.gz \\
        step1.vcf.gz
    
    # Apply depth filter
    bcftools view \\
        -i "FORMAT/DP>=5" \\
        -O z \\
        -o step3_depth.vcf.gz \\
        step2.vcf.gz
    
    # Apply quality filter
    if zcat step3_depth.vcf.gz | grep -E "^##FORMAT=<ID=GQ|^##INFO=<ID=GQ" >/dev/null; then
        bcftools view \\
            -i "FORMAT/GQ>=20 || INFO/GQ>=20 || QUAL>=20" \\
            -O v \\
            step3_depth.vcf.gz | \\
        bcftools sort -O z -o ${sample_id}.filtered_variants.vcf.gz
    else
        bcftools view \\
            -i "QUAL>=20" \\
            -O v \\
            step3_depth.vcf.gz | \\
        bcftools sort -O z -o ${sample_id}.filtered_variants.vcf.gz
    fi
    
    # Index the filtered VCF
    tabix -p vcf ${sample_id}.filtered_variants.vcf.gz
    
    rm step1.vcf.gz step2.vcf.gz step3_depth.vcf.gz
    """
}

process annotate {

    input:
    tuple val(sample_id), path(vcf_file)
    
    output:
    path "${sample_id}.annotated_variants.vcf.gz"

    script:
    def annotation_table = "${baseDir}/data/enhanced_annotation_table_4.tsv.gz"
    def annotation_header = "${baseDir}/data/enhanced_annotations_header_4.txt"  
    """
    bcftools annotate \\
        -a ${annotation_table} \\
        -h ${annotation_header} \\
        -c CHROM,POS,REF,ALT,INFO/GENE,INFO/DRUG,INFO/EFFECT,INFO/WHO_CLASSIFICATION,INFO/VARIANT_ID,INFO/GENOME_POSITION \\
        -O z \\
        -o ${sample_id}.annotated_variants.vcf.gz \\
        ${vcf_file}
    """
}

workflow VCF_PROCESSING {
    take:
    vcf_ch  

    main:
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
    normalized = VCF_NORMALIZE(vcf_ch, reference_ch.first())
    filtered = filter_variants(normalized)
    annotated = annotate(filtered)

    emit:
    annotated = annotated           
    filtered = filtered           

}