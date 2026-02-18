#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info """
    Mycobacterium tuberculosis Mutation Analysis Pipeline (v${params.version})
    Developed by SPHERES-OUCRU ID
    Documentation: https://tb-pipeline-docs.readthedocs.io/
"""

include { ILLUMINA }              from './workflows/illumina.nf'
include { NANOPORE }              from './workflows/nanopore.nf'
include { VCF_PROCESSING }        from './workflows/vcf.nf'  
include { LINEAGE }               from './workflows/lineage.nf'
include { GENERATE_REPORT }       from './workflows/report.nf'
include { GENERATE_SAMPLE_REPORTS } from './workflows/report.nf'
include { FHIR }                  from './workflows/fhir.nf'
include { VALIDATE }              from './workflows/validate_fhir.nf'
include { MERGE_CLINICAL_DATA }   from './workflows/merge_clinical_data.nf'
include { UPLOAD_FHIR }           from './workflows/upload_fhir.nf'
include { VERSIONS }              from './workflows/utils.nf'
include { DEEPLEX }               from './workflows/deeplex.nf'

workflow {
    illumina_reads_ch = Channel
        .fromFilePairs("${params.reads_dir}/*_{1,2}_illumina.fastq.gz")
        .map { id, files -> tuple(id, files) }  

    nanopore_reads_ch = Channel
        .fromPath("${params.reads_dir}/*_ont.fastq.gz", checkIfExists: false)
        .map { file -> tuple(file.baseName.replaceFirst(/_ont$/, ''), file) }

    vcf_ch = Channel
        .fromPath("${params.vcf_dir}/*.vcf{,.gz}", checkIfExists: false)
        .map { file -> tuple(file.baseName.replaceFirst(/\.vcf(\.gz)?$/, ''), file) }

    deeplex_ch = Channel
        .fromPath("${params.deeplex_dir}/*.xlsx", checkIfExists: false)
        .filter { file -> !file.name.startsWith('~$') }

    illumina_out = ILLUMINA(illumina_reads_ch)
    nanopore_out = NANOPORE(nanopore_reads_ch)
    vcf_out = VCF_PROCESSING(vcf_ch) 

    // all_consensus = illumina_out.consensus
    //    .mix(nanopore_out.consensus)
    //    .collect()

    all_filtered = illumina_out.filtered
        .mix(nanopore_out.filtered)
        .mix(vcf_out.filtered) 

    lineage_out = LINEAGE(all_filtered)

    all_qc = illumina_out.qc_report
        .mix(nanopore_out.qc_report)
        .map { it -> it[1] }
        .flatten()
        .filter { it.toString().endsWith('_fastqc.zip') }
        .collect()

    all_annotated = illumina_out.annotated
        .mix(nanopore_out.annotated)
        .mix(vcf_out.annotated)  
        .map { it -> it instanceof List ? it[1] : it }
        
    GENERATE_REPORT(all_qc)
    
    lineage_files = lineage_out.lineage_results
        .map { sample_id, file_path -> file_path } 
        .collect()
    
    sample_reports = GENERATE_SAMPLE_REPORTS(all_annotated, lineage_files)
    
    fhir_out = FHIR(all_annotated, lineage_files) 
    
    clinical_metadata_ch = Channel.fromPath(params.clinical_metadata, checkIfExists: false)
        .first() 
    
    deeplex_clinical_ch = Channel.fromPath(params.clinical_metadata_deeplex, checkIfExists: false)
        .first()

    deeplex_out = DEEPLEX(deeplex_ch, deeplex_clinical_ch)

    merged_clinical_out = MERGE_CLINICAL_DATA(
        fhir_out.fhir_output, 
        clinical_metadata_ch
    )
    
    validation_out = VALIDATE(merged_clinical_out.merged_fhir)
    
    // Optional: Upload validated FHIR
    //upload_out = UPLOAD_FHIR(validation_out.validated_fhir)

    VERSIONS()
}