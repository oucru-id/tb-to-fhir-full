#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process fastqc {

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), emit: qc_report  

    script:
    """
    fastqc ${reads} --outdir .
    """
}

process chopper {
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fastq.gz"), emit: trimmed_reads

    script:
    """
    zcat ${reads} | chopper -q ${params.nanopore_min_q} -l ${params.nanopore_min_l} | gzip > ${sample_id}_trimmed.fastq.gz
    """
}

process minimap2 {

    input:
    tuple val(sample_id), path(reads), path(reference)

    output:
    tuple val(sample_id), path("aligned.bam")

    script:
    """

    minimap2 -x map-ont --secondary=no -L --MD \\
        -A 2 -B 4 -O 4,24 -E 2,1 \\
        -t ${task.cpus} -a ${reference} ${reads} \\
    | samtools sort -o aligned.bam
    samtools index aligned.bam
    """
}

process medaka {

    input:
    tuple val(sample_id), path("aligned.bam"), path(reference)

    output:
    tuple val(sample_id), path("variants.vcf.gz")

    script:

    def model = params.medaka_model ?: 'r941_e81_sup_variant_g514'
    def batch_size = params.medaka_batch_size ?: 100
    """
    if [ -L "${reference}" ]; then
        actual_reference=\$(readlink -f "${reference}")
        cp "\${actual_reference}" reference.fasta
    else
        cp "${reference}" reference.fasta
    fi
    reference_path="reference.fasta"

    if [ ! -s aligned.bam ]; then
        ls -lh aligned.bam
        exit 1
    fi

    samtools index aligned.bam
    samtools faidx \$reference_path

    export CUDA_VISIBLE_DEVICES=""

    medaka inference aligned.bam consensus_probs.hdf \\
        --model ${model} \\
        --batch_size ${batch_size} \\
        --threads ${task.cpus} \\
        || { echo "Failed to run medaka inference."; exit 1; }

    medaka vcf consensus_probs.hdf \$reference_path medaka.vcf \\
        || { echo "Failed to create variants."; exit 1; }

    bcftools sort medaka.vcf -o medaka.sorted.vcf \\
        || { echo "Failed to sort variants."; exit 1; }

    medaka tools annotate medaka.sorted.vcf \$reference_path aligned.bam medaka.annotated.vcf \\
        || { echo "Failed to annotate variants."; exit 1; }

    if [ ! -s medaka.annotated.vcf ]; then
        echo "medaka produced no annotated VCF"
        ls -la
        exit 1
    fi

    bgzip -c medaka.annotated.vcf > variants.vcf.gz
    """

    stub:
    """
    touch variants.vcf.gz
    """
}

process annotate {

    input:
    tuple val(sample_id), path("variants.vcf.gz")
    
    output:
    path "${sample_id}.annotated_variants.vcf.gz"

    script:
    def annotation_table = "${params.annotation_table}"
    def annotation_header = "${params.annotation_header}"
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
        -i 'DP>=0 && (GT="0" || GT="1" || GT="0/1")' \
        -O z \
        -o step2.vcf.gz \
        step1.vcf.gz
    
    # Apply depth filter
    bcftools view \
        -i "DP>=5" \
        -O z \
        -o step3_depth.vcf.gz \
        step2.vcf.gz
    
    # Apply quality filter
    bcftools view \
        -i "GQ>=20" \
        -O v \
        step3_depth.vcf.gz | \
    bcftools sort -O z -o ${sample_id}.filtered_variants.vcf.gz
    
    # Index the filtered VCF
    tabix -p vcf ${sample_id}.filtered_variants.vcf.gz
    
    # Clean up intermediate files
    rm step1.vcf.gz step2.vcf.gz step3_depth.vcf.gz
    
    """
}

workflow NANOPORE {
    take:
    reads

    main:
    qc = fastqc(reads)
    trimmed = chopper(reads)
    aligned = minimap2(trimmed.trimmed_reads.map { it -> tuple(it[0], it[1], file(params.reference)) })

    variants = medaka(aligned.map { it -> tuple(it[0], it[1], file(params.reference)) })
    filtered = filter_variants(variants)
    annotated = annotate(filtered)

    emit:
    qc_report = qc.qc_report
    annotated = annotated
    filtered = filtered
    bam = aligned
}
