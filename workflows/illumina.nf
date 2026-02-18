#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process fastqc {

    input:
    tuple val(sample_id), path(reads)  

    output:
    tuple val(sample_id), path("*_fastqc.zip"), emit: qc_report

    script:
    """
    fastqc ${reads.join(' ')} --outdir .
    """
}

process trimmomatic {

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(paired_out), emit: trimmed_reads

    script:
    def r1 = reads[0]
    def r2 = reads[1]
    paired_out = ["${sample_id}_1.trimmed.fastq.gz", "${sample_id}_2.trimmed.fastq.gz"]
    
    """
    java -jar /usr/share/java/trimmomatic.jar PE -threads ${task.cpus} -phred33 \\
        ${r1} ${r2} \\
        ${sample_id}_1.trimmed.fastq.gz ${sample_id}_1.unpaired.fastq.gz \\
        ${sample_id}_2.trimmed.fastq.gz ${sample_id}_2.unpaired.fastq.gz \\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
    """
}

process bwa_mem2 {

    input:
    tuple val(sample_id), path(reads), path(reference)

    output:
    tuple val(sample_id), path("aligned.bam")

    script:
    """
    if [ ! -f ${reference}.bwt.2bit.64 ]; then
        bwa-mem2 index ${reference}
    fi

    bwa-mem2 mem -t ${task.cpus} \
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:lib1" \
        ${reference} \
        ${reads[0]} ${reads[1]} | \
    samtools sort -o aligned.bam

    samtools index aligned.bam
    """
}

process gatk {

    input:
    tuple val(sample_id), path("aligned.bam"), path(reference)

    output:
    tuple val(sample_id), path("variants.vcf.gz")

    script:
    def gatk_path = "${baseDir}/tools/gatk-4.4.0.0/gatk"
    def ref_prefix = reference.getBaseName()
    """
    if [ ! -f ${reference}.fai ]; then
        samtools faidx ${reference}
    fi
    if [ ! -f ${ref_prefix}.dict ]; then
        ${gatk_path} CreateSequenceDictionary -R ${reference}
    fi

    cp aligned.bam clean.bam
    samtools index clean.bam

    ${gatk_path} HaplotypeCaller \
        -R ${reference} \
        -I clean.bam \
        -O variants.vcf.gz \
        --create-output-variant-index true \
        --native-pair-hmm-threads ${task.cpus} \
        --min-base-quality-score 10 \
        --standard-min-confidence-threshold-for-calling 10.0 \
        --annotation MappingQuality \
        --annotation QualByDepth \
        --annotation FisherStrand \
        --output-mode EMIT_VARIANTS_ONLY \
        --java-options "-Xmx4g" \
        --dont-use-soft-clipped-bases \
        --assembly-region-padding 100 \
        --max-reads-per-alignment-start 0 \
        --min-pruning 1

    if [ ! -f variants.vcf.gz.tbi ]; then
        tabix -p vcf variants.vcf.gz
    fi
    """
}

process annotate {

    input:
    tuple val(sample_id), path("variants.vcf.gz")
    
    output:
    path "${sample_id}.annotated_variants.vcf.gz"

    script:
    def annotation_table = "${baseDir}/data/enhanced_annotation_table_4.tsv.gz"
    def annotation_header = "${baseDir}/data/enhanced_annotations_header_4.txt"  
    """
    bcftools annotate \
        -a ${annotation_table} \
        -h ${annotation_header} \
        -c CHROM,POS,REF,ALT,INFO/GENE,INFO/DRUG,INFO/EFFECT,INFO/WHO_CLASSIFICATION,INFO/VARIANT_ID,INFO/GENOME_POSITION \
        -O z \
        -o ${sample_id}.annotated_variants.vcf.gz \
        variants.vcf.gz
    """
}

process filter_variants {

    input:
    tuple val(sample_id), path("variants.vcf.gz")
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered_variants.vcf.gz")

    script:
    def repetitive_regions = "${baseDir}/data/repetitive_regions.bed"
    """
    # Exclude repetitive regions (pe/ppe genes, insertion sequences)
    if [ -f ${repetitive_regions} ]; then
        bcftools view -T ^${repetitive_regions} variants.vcf.gz -O z -o step1.vcf.gz
    else
        cp variants.vcf.gz step1.vcf.gz
    fi
    
    # Filtering 
    bcftools view \
        -v snps,indels \
        -i 'FORMAT/DP>=0 && (GT="0/0" || GT="1/1" || GT="0/1")' \
        -O z \
        -o step2.vcf.gz \
        step1.vcf.gz
    
    bcftools view \
        -i "FORMAT/DP>=5" \
        -O z \
        -o step3_depth.vcf.gz \
        step2.vcf.gz
    
    # Apply quality-based filter
    if zcat step3_depth.vcf.gz | grep "^##FORMAT=<ID=GQ" >/dev/null; then
        bcftools view \
            -i "FORMAT/GQ>=20" \
            -O v \
            step3_depth.vcf.gz | \
        bcftools sort -O z -o ${sample_id}.filtered_variants.vcf.gz
    else
        bcftools sort -O z -o ${sample_id}.filtered_variants.vcf.gz step3_depth.vcf.gz
    fi
    
    # Index the filtered VCF
    tabix -p vcf ${sample_id}.filtered_variants.vcf.gz
    
    # Clean up intermediate files
    rm step1.vcf.gz step2.vcf.gz step3_depth.vcf.gz
    
    """
}

workflow ILLUMINA {
    take:
    reads

    main:
    qc = fastqc(reads)
    
    trimmed = trimmomatic(reads)
    
    aligned = bwa_mem2(trimmed.trimmed_reads.map { id, reads -> tuple(id, reads, file(params.reference)) })
        
    variants = gatk(aligned.map { id, bam -> tuple(id, bam, file(params.reference)) })
    filtered = filter_variants(variants)
    annotated = annotate(filtered)  

    emit:
    qc_report = qc.qc_report
    annotated = annotated
    filtered = filtered
}