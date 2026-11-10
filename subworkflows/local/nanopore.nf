include { FASTQC         } from '../../modules/local/fastqc'
include { CHOPPER        } from '../../modules/local/chopper'
include { MINIMAP2_ALIGN } from '../../modules/local/minimap2_align'
include { MEDAKA_VARIANT } from '../../modules/local/medaka_variant'

workflow NANOPORE {
    take:
    reads      
    fasta

    main:
    FASTQC(reads)
    CHOPPER(reads)
    MINIMAP2_ALIGN(CHOPPER.out.reads, fasta)
    MEDAKA_VARIANT(MINIMAP2_ALIGN.out.bam, fasta)

    emit:
    vcf      = MEDAKA_VARIANT.out.vcf
    bam      = MINIMAP2_ALIGN.out.bam
    qc       = FASTQC.out.zip
    versions = FASTQC.out.versions.mix(CHOPPER.out.versions, MINIMAP2_ALIGN.out.versions,
                                       MEDAKA_VARIANT.out.versions)
}
