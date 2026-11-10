include { FASTQC                } from '../../modules/local/fastqc'
include { TRIMMOMATIC           } from '../../modules/local/trimmomatic'
include { BWAMEM2_MEM           } from '../../modules/local/bwamem2_mem'
include { GATK4_HAPLOTYPECALLER } from '../../modules/local/gatk4_haplotypecaller'

workflow ILLUMINA {
    take:
    reads 
    fasta
    fai
    dict
    bwa_index

    main:
    FASTQC(reads)
    TRIMMOMATIC(reads)
    BWAMEM2_MEM(TRIMMOMATIC.out.reads, bwa_index)
    GATK4_HAPLOTYPECALLER(BWAMEM2_MEM.out.bam, fasta, fai, dict)

    emit:
    vcf      = GATK4_HAPLOTYPECALLER.out.vcf
    bam      = BWAMEM2_MEM.out.bam
    qc       = FASTQC.out.zip
    versions = FASTQC.out.versions.mix(TRIMMOMATIC.out.versions, BWAMEM2_MEM.out.versions,
                                       GATK4_HAPLOTYPECALLER.out.versions)
}
