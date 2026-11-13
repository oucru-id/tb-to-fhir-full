include { EXTRACT_BARCODE_SNPS } from '../../modules/local/extract_barcode_snps'
include { CLASSIFY_LINEAGE     } from '../../modules/local/classify_lineage'

workflow LINEAGE {
    take:
    filtered 
    barcode_bed

    main:
    EXTRACT_BARCODE_SNPS(filtered, barcode_bed)
    CLASSIFY_LINEAGE(EXTRACT_BARCODE_SNPS.out.vcf, barcode_bed)

    emit:
    json     = CLASSIFY_LINEAGE.out.json
    versions = EXTRACT_BARCODE_SNPS.out.versions.mix(CLASSIFY_LINEAGE.out.versions)
}
