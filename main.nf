#!/usr/bin/env nextflow

include { paramsSummaryLog    } from 'plugin/nf-schema'
include { validateTypedParams } from './subworkflows/local/utils_params'
include { TBFHIR              } from './workflows/tbfhir'

workflow {
    main:
    log.info """\
        TBtoFHIR ${workflow.manifest.version}: Mycobacterium tuberculosis genomics to HL7 FHIR
        ${workflow.manifest.homePage}
        """.stripIndent()
    if (params.results_dir) {
        error "--results_dir was renamed to --outdir in 2.0.0"
    }
    if (workflow.resume && nextflow.version.toString() in ['26.04.0', '26.04.1']) {
        log.warn "Nextflow ${nextflow.version} does not re-run a task when a parameter it uses has changed, " +
                 "so -resume can reuse results made with other settings. " +
                 "Update Nextflow, or run without -resume after changing parameters."
    }
    validateTypedParams()
    log.info paramsSummaryLog(workflow)
    TBFHIR()
}
