
include { SAMTOOLS_FAIDX                 } from '../../modules/local/samtools_faidx'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/local/gatk4_createsequencedictionary'
include { BWAMEM2_INDEX                  } from '../../modules/local/bwamem2_index'

workflow PREPARE_REFERENCE {
    take:
    fasta          
    has_illumina   

    main:
    ch_versions = channel.empty()
    SAMTOOLS_FAIDX(fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    ch_dict = channel.empty()
    ch_bwa  = channel.empty()
    if (has_illumina) {
        GATK4_CREATESEQUENCEDICTIONARY(fasta)
        BWAMEM2_INDEX(fasta)
        ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        ch_bwa  = BWAMEM2_INDEX.out.index
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions, BWAMEM2_INDEX.out.versions)
    }

    emit:
    fai       = SAMTOOLS_FAIDX.out.fai    
    dict      = ch_dict
    bwa_index = ch_bwa
    versions  = ch_versions
}
