
include { FHIR_GENOMIC                           } from '../../modules/local/fhir_genomic'
include { FHIR_CLINICAL                          } from '../../modules/local/fhir_clinical'
include { FHIR_LINK                              } from '../../modules/local/fhir_link'
include { FHIR_VALIDATE; FHIR_VALIDATION_SUMMARY } from '../../modules/local/fhir_validate'
include { FHIR_UPLOAD; FHIR_UPLOAD_WITH_API_KEY } from '../../modules/local/fhir_upload'

workflow FHIR {
    take:
    annotated    
    lineage      
    coverage     
    gene_map
    no_file

    main:
    ch_versions = channel.empty()
    ch_genomic_in = annotated
        .map { meta, vcf -> [meta.id, meta, vcf] }
        .join(lineage.map { meta, json -> [meta.id, json] }, failOnMismatch: true, failOnDuplicate: true)
        .join(coverage.map { meta, json -> [meta.id, json] }, failOnMismatch: true, failOnDuplicate: true)
        .map { _id, meta, vcf, lineage_json, coverage_json -> [meta, vcf, lineage_json, coverage_json] }

    FHIR_GENOMIC(ch_genomic_in, gene_map)
    ch_versions = ch_versions.mix(FHIR_GENOMIC.out.versions)
    ch_genomic = FHIR_GENOMIC.out.json.map { _meta, json -> json }

    ch_bundles = channel.empty()
    if (params.patient_metadata) {
        ch_ids = FHIR_GENOMIC.out.json
            .map { meta, _json -> meta.id }
            .collectFile(name: 'sample_ids.txt', newLine: true, sort: true)
        FHIR_CLINICAL(
            file(params.patient_metadata, checkIfExists: true),
            params.organization_metadata ? file(params.organization_metadata, checkIfExists: true) : no_file,
            params.practitioner_metadata ? file(params.practitioner_metadata, checkIfExists: true) : no_file,
            ch_ids)
        FHIR_LINK(
            ch_genomic.collect(),
            FHIR_CLINICAL.out.dir,
            params.external_refs ? file(params.external_refs, checkIfExists: true) : no_file)
        ch_bundles = FHIR_LINK.out.bundles.flatten()
        ch_versions = ch_versions.mix(FHIR_CLINICAL.out.versions, FHIR_LINK.out.versions)
    }

    ch_uploadable = ch_bundles
    if (!params.skip_fhir_validation) {
        FHIR_VALIDATE(
            ch_genomic.mix(ch_bundles).toSortedList().flatMap { files -> files.collate(params.fhir_validation_batch_size) },
            file(params.fhir_validator),
            file("${projectDir}/assets/fhir", checkIfExists: true))
        FHIR_VALIDATION_SUMMARY(FHIR_VALIDATE.out.result)
        ch_versions = ch_versions.mix(FHIR_VALIDATE.out.versions.first(), FHIR_VALIDATION_SUMMARY.out.versions.first())
        ch_uploadable = FHIR_VALIDATION_SUMMARY.out.passed.flatten().filter { f -> f.name.endsWith('.bundle.fhir.json') }
    }

    if (params.fhir_server_url) {
        if (params.fhir_upload_api_key) {
            FHIR_UPLOAD_WITH_API_KEY(ch_uploadable)
            ch_versions = ch_versions.mix(FHIR_UPLOAD_WITH_API_KEY.out.versions.first())
        } else {
            FHIR_UPLOAD(ch_uploadable)
            ch_versions = ch_versions.mix(FHIR_UPLOAD.out.versions.first())
        }
    }

    emit:
    genomic  = FHIR_GENOMIC.out.json
    bundles  = ch_bundles
    versions = ch_versions
}
