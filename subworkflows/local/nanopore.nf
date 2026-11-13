include { NANOPLOT       } from '../../modules/local/nanoplot'
include { CHOPPER        } from '../../modules/local/chopper'
include { MINIMAP2_ALIGN } from '../../modules/local/minimap2_align'
include { MEDAKA_VARIANT } from '../../modules/local/medaka_variant'

workflow NANOPORE {
    take:
    reads      
    fasta

    main:
    NANOPLOT(reads)
    CHOPPER(reads)
    MINIMAP2_ALIGN(CHOPPER.out.reads, fasta)
    MEDAKA_VARIANT(MINIMAP2_ALIGN.out.bam, fasta)

    emit:
    vcf      = MEDAKA_VARIANT.out.vcf
    bam      = MINIMAP2_ALIGN.out.bam
    qc       = NANOPLOT.out.stats
    versions = NANOPLOT.out.versions.mix(CHOPPER.out.versions, MINIMAP2_ALIGN.out.versions,
                                         MEDAKA_VARIANT.out.versions)
}
