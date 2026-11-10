#!/usr/bin/env nextflow

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { TBFHIR } from './workflows/tbfhir'

workflow {
    main:
    log.info """\
        TBtoFHIR ${workflow.manifest.version}: Mycobacterium tuberculosis genomics to HL7 FHIR
        ${workflow.manifest.homePage}
        """.stripIndent()
    if (params.results_dir) {
        error "--results_dir was renamed to --outdir in 2.0.0"
    }
    validateParameters()
    log.info paramsSummaryLog(workflow)
    TBFHIR()
}
