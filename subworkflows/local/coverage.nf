
include { MAKE_WHO_TARGETS                          } from '../../modules/local/make_who_targets'
include { MOSDEPTH                                  } from '../../modules/local/mosdepth'
include { SUMMARIZE_COVERAGE; COVERAGE_NOT_ASSESSED } from '../../modules/local/summarize_coverage'

workflow COVERAGE {
    take:
    bams                 
    vcf_only_meta        
    annotation_table
    repetitive_regions
    fai

    main:
    MAKE_WHO_TARGETS(annotation_table, repetitive_regions, fai)
    targets = MAKE_WHO_TARGETS.out.bed
    MOSDEPTH(bams, targets)
    SUMMARIZE_COVERAGE(MOSDEPTH.out.beds, targets)
    COVERAGE_NOT_ASSESSED(vcf_only_meta)

    emit:
    json     = SUMMARIZE_COVERAGE.out.json.mix(COVERAGE_NOT_ASSESSED.out.json)
    versions = MAKE_WHO_TARGETS.out.versions.mix(MOSDEPTH.out.versions, SUMMARIZE_COVERAGE.out.versions,
                                                 COVERAGE_NOT_ASSESSED.out.versions)
}
