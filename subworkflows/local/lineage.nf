include { BARCODE_PILEUP   } from '../../modules/local/barcode_pileup'
include { CLASSIFY_LINEAGE } from '../../modules/local/classify_lineage'

workflow LINEAGE {
    take:
    filtered
    bams
    fasta
    fai
    barcode_bed

    main:
    BARCODE_PILEUP(bams, fasta, fai, barcode_bed)
    ch_classify = filtered
        .map { meta, vcf, _tbi -> [meta.id, meta, vcf] }
        .join(BARCODE_PILEUP.out.tsv.map { meta, tsv -> [meta.id, tsv] }, remainder: true)
        .map { _id, meta, vcf, tsv -> [meta, vcf, tsv ?: []] }
    CLASSIFY_LINEAGE(ch_classify, fasta, barcode_bed)

    emit:
    json     = CLASSIFY_LINEAGE.out.json
    versions = BARCODE_PILEUP.out.versions.mix(CLASSIFY_LINEAGE.out.versions)
}
