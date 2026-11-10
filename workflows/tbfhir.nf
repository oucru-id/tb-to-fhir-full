
include { samplesheetToList } from 'plugin/nf-schema'
include { PREPARE_REFERENCE } from '../subworkflows/local/prepare_reference'
include { ILLUMINA          } from '../subworkflows/local/illumina'
include { NANOPORE          } from '../subworkflows/local/nanopore'
include { LINEAGE           } from '../subworkflows/local/lineage'
include { COVERAGE          } from '../subworkflows/local/coverage'
include { FHIR              } from '../subworkflows/local/fhir'
include { VCF_NORMALIZE     } from '../modules/local/vcf_normalize'
include { FILTER_VARIANTS   } from '../modules/local/filter_variants'
include { ANNOTATE_VARIANTS } from '../modules/local/annotate_variants'
include { MULTIQC           } from '../modules/local/multiqc'
include { SAMPLE_REPORT     } from '../modules/local/sample_report'

def fhirId(String id) {
    return id.replaceAll(/[^A-Za-z0-9\-]/, '-').replaceAll(/-+/, '-').replaceAll(/^-+|-+$/, '')
}

def present(value) {
    return value != null && value != [] && value.toString() != ''
}

def samplesFromSheet(String sheet) {
    def samples = []
    def rows = samplesheetToList(sheet, "${projectDir}/assets/schema_input.json")
    rows.eachWithIndex { row, index ->
        def meta = row[0]
        def fastq_1 = row[1]
        def fastq_2 = row[2]
        def vcf = row[3]
        def line = index + 2
        if (meta.platform == 'illumina') {
            if (!present(fastq_1) || !present(fastq_2)) {
                error "samplesheet line ${line} (${meta.id}): illumina needs fastq_1 and fastq_2"
            }
            samples << [meta: [id: meta.id, platform: 'illumina'], files: [fastq_1, fastq_2]]
        } else if (meta.platform == 'ont') {
            if (!present(fastq_1)) {
                error "samplesheet line ${line} (${meta.id}): ont needs fastq_1"
            }
            samples << [meta: [id: meta.id, platform: 'ont'], files: fastq_1]
        } else {
            if (!present(vcf)) {
                error "samplesheet line ${line} (${meta.id}): vcf needs a vcf column"
            }
            samples << [meta: [id: meta.id, platform: 'vcf'], files: vcf]
        }
    }
    return samples
}


def samplesFromDirectories() {
    def samples = []
    if (params.reads_dir) {
        def pairs = [:]
        files("${params.reads_dir}/*_{1,2}_illumina.fastq.gz").each { f ->
            def id = f.name.replaceFirst(/_[12]_illumina\.fastq\.gz$/, '')
            if (!pairs.containsKey(id)) {
                pairs[id] = []
            }
            pairs[id] << f
        }
        pairs.keySet().sort().each { id ->
            def reads = pairs[id].sort { f -> f.name }
            if (reads.size() != 2) {
                error "Illumina sample ${id}: expected _1_illumina and _2_illumina files, found ${reads*.name}"
            }
            samples << [meta: [id: id, platform: 'illumina'], files: reads]
        }
        files("${params.reads_dir}/*_ont.fastq.gz").sort { f -> f.name }.each { f ->
            samples << [meta: [id: f.name.replaceFirst(/_ont\.fastq\.gz$/, ''), platform: 'ont'], files: f]
        }
    }
    if (params.vcf_dir) {
        def vcfs = files("${params.vcf_dir}/*.vcf") + files("${params.vcf_dir}/*.vcf.gz")
        vcfs.sort { f -> f.name }.each { f ->
            samples << [meta: [id: f.name.replaceFirst(/\.vcf(\.gz)?$/, ''), platform: 'vcf'], files: f]
        }
    }
    return samples
}

def validateSamples(samples) {
    if (!samples) {
        error "No samples found. Give --input samplesheet.csv, or --reads_dir / --vcf_dir " +
              "(files named <id>_1_illumina.fastq.gz, <id>_ont.fastq.gz or <id>.vcf.gz)."
    }
    def seen = [:]
    def seenFhir = [:]
    samples.each { s ->
        def id = s.meta.id
        if (!(id ==~ /^[A-Za-z0-9][A-Za-z0-9._-]{0,39}$/)) {
            error "Sample id '${id}': use 1-40 letters, digits, '.', '_' or '-', starting with a letter or digit"
        }
        if (seen.containsKey(id)) {
            error "Sample id '${id}' appears more than once"
        }
        def fid = fhirId(id)
        if (seenFhir.containsKey(fid)) {
            error "Sample ids '${seenFhir[fid]}' and '${id}' become the same FHIR id '${fid}'; rename one"
        }
        if (id ==~ /^\d{16}$/) {
            log.warn "Sample id '${id}' looks like a 16-digit national ID (NIK). Sample ids appear in " +
                     "the PHI-free genomic output; use a laboratory accession instead."
        }
        seen[id] = true
        seenFhir[fid] = id
    }
    if (samples.any { s -> s.meta.platform == 'ont' } && !params.medaka_model) {
        error "ONT samples need --medaka_model (e.g. r1041_e82_400bps_sup_variant_v5.0.0); " +
              "list the available models with `medaka tools list_models`."
    }
    if (params.fhir_server_url && !params.patient_metadata) {
        error "--fhir_server_url uploads linked bundles, which need --patient_metadata"
    }
    if (params.repetitive_regions && file(params.repetitive_regions).name == 'NO_FILE') {
        log.warn "Repetitive-region masking is OFF (repetitive_regions = NO_FILE)"
    }
}


workflow TBFHIR {
    main:
    def samples = params.input ? samplesFromSheet(params.input) : samplesFromDirectories()
    validateSamples(samples)
    def counts = samples.countBy { s -> s.meta.platform }
    log.info "Samples: ${samples.size()} (illumina: ${counts.illumina ?: 0}, ont: ${counts.ont ?: 0}, vcf: ${counts.vcf ?: 0})"

    def no_file = file("${projectDir}/assets/NO_FILE", checkIfExists: true)
    def fasta = file(params.reference, checkIfExists: true)
    def annotation_table = file(params.annotation_table, checkIfExists: true)
    def annotation_index = file("${params.annotation_table}.tbi", checkIfExists: true)
    def annotation_header = file(params.annotation_header, checkIfExists: true)
    def barcode_bed = file(params.lineage_barcode, checkIfExists: true)
    def repetitive = params.repetitive_regions ? file(params.repetitive_regions, checkIfExists: true) : no_file
    def gene_map = params.gene_map ? file(params.gene_map, checkIfExists: true) : no_file

    ch_samples  = channel.fromList(samples).map { s -> [s.meta, s.files] }
    ch_illumina = ch_samples.filter { meta, _f -> meta.platform == 'illumina' }
    ch_ont      = ch_samples.filter { meta, _f -> meta.platform == 'ont' }
    ch_vcf      = ch_samples.filter { meta, _f -> meta.platform == 'vcf' }
    ch_versions = channel.empty()

    PREPARE_REFERENCE(fasta, (counts.illumina ?: 0) > 0)
    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    ILLUMINA(ch_illumina, fasta, PREPARE_REFERENCE.out.fai, PREPARE_REFERENCE.out.dict,
             PREPARE_REFERENCE.out.bwa_index)
    NANOPORE(ch_ont, fasta)
    VCF_NORMALIZE(ch_vcf, fasta, PREPARE_REFERENCE.out.fai)
    ch_versions = ch_versions.mix(ILLUMINA.out.versions, NANOPORE.out.versions, VCF_NORMALIZE.out.versions)

    FILTER_VARIANTS(ILLUMINA.out.vcf.mix(NANOPORE.out.vcf, VCF_NORMALIZE.out.vcf), repetitive)
    ANNOTATE_VARIANTS(FILTER_VARIANTS.out.vcf, annotation_table, annotation_index, annotation_header)
    ch_versions = ch_versions.mix(FILTER_VARIANTS.out.versions, ANNOTATE_VARIANTS.out.versions)

    LINEAGE(FILTER_VARIANTS.out.vcf, barcode_bed)
    COVERAGE(ILLUMINA.out.bam.mix(NANOPORE.out.bam), ch_vcf.map { meta, _f -> meta },
             annotation_table, repetitive, PREPARE_REFERENCE.out.fai)
    ch_versions = ch_versions.mix(LINEAGE.out.versions, COVERAGE.out.versions)

    ch_qc = ILLUMINA.out.qc.mix(NANOPORE.out.qc).flatMap { _meta, zips -> [zips].flatten() }
    MULTIQC(ch_qc.collect())
    ch_versions = ch_versions.mix(MULTIQC.out.versions)

    SAMPLE_REPORT(
        ANNOTATE_VARIANTS.out.vcf
            .map { meta, vcf -> [meta.id, meta, vcf] }
            .join(LINEAGE.out.json.map { meta, json -> [meta.id, json] }, failOnMismatch: true)
            .map { _id, meta, vcf, json -> [meta, vcf, json] })
    ch_versions = ch_versions.mix(SAMPLE_REPORT.out.versions)

    FHIR(ANNOTATE_VARIANTS.out.vcf, LINEAGE.out.json, COVERAGE.out.json, gene_map, no_file)
    ch_versions = ch_versions.mix(FHIR.out.versions)

    def versions_header = ("\"Workflow\":\n    TBtoFHIR: ${workflow.manifest.version}\n" +
                           "    Nextflow: ${nextflow.version}\n" +
                           "    fhir_profile: ${params.fhir_profile}\n" +
                           "    reference: ${fasta.name}\n    annotation_table: ${annotation_table.name}\n" +
                           "    lineage_barcode: ${barcode_bed.name}\n").toString()
    ch_versions
        .unique()
        .collectFile(name: 'versions_collected.yml', sort: true, newLine: true)
        .map { collected -> versions_header + collected.text }
        .collectFile(name: 'software_versions.yml', storeDir: "${params.outdir}/pipeline_info")

    emit:
    genomic = FHIR.out.genomic
    bundles = FHIR.out.bundles
}
