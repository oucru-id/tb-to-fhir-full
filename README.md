# TB Genomics Mutation Analysis to FHIR Genomics Pipeline (TBtoFHIR)

A Nextflow pipeline for *Mycobacterium tuberculosis* genomic analysis from raw sequencing data or
VCFs. It produces WHO-catalogue drug-resistance calls, lineage and coverage QC, and reports them as
HL7 FHIR R4 resources conforming to the
[Genomics Reporting IG 3.0.0](https://hl7.org/fhir/uv/genomics-reporting/STU3/).

## Key features

- **Multi-platform input**: Illumina paired-end reads, Oxford Nanopore reads, or VCFs.
- **Per-drug WHO grading**
- **Coverage-gated susceptibility**: a drug is reported susceptible only when every one of its
  resistance loci was sequenced to the required depth and breadth. Otherwise it is
  **Indeterminate**.
- **Lineage classification** from barcode SNPs, with confidence and SNP support.
- **QC**: FastQC (Illumina) or NanoPlot (ONT), MultiQC, and per-locus coverage over the WHO
  resistance targets.
- **Two FHIR profiles**:
  - `generic` (default): plain HL7 Genomics Reporting.
  - `satusehat`: adds what the Indonesian Ministry of Health SATUSEHAT platform expects (Kemkes
    profiles, NIK identifiers, Kemkes codings).

## Requirements

You need [Nextflow](https://www.nextflow.io/) ≥ 25.04 and **one** of the following. On 26.04, use
26.04.2 or later.

| Profile | Provides the tools via |
|---|---|
| `-profile docker` (or `podman`) | pinned [BioContainers](https://biocontainers.pro/) images |
| `-profile singularity` (or `apptainer`) | the same images, e.g. on HPC |
| `-profile conda` (or `mamba`) | pinned bioconda / conda-forge packages |
| none | tools already on `PATH` (see below) |

- FastQC 0.12.1, Trimmomatic 0.39
- bwa-mem2 2.2.1, minimap2 2.26, samtools 1.19.2
- GATK 4.4.0.0, medaka 2.1.0, chopper 0.12.0
- bcftools 1.19, mosdepth 0.3.6, MultiQC 1.30
- Python 3.12, Java 21
- NanoPlot 1.48.0, the ONT read QC, is new in 2.0.

Without a profile you need on `PATH`:
- for all runs: `bcftools`, `samtools`, `tabix`/`bgzip`, Python 3 (standard library only) and
  Java ≥ 17
- for Illumina input: `bwa-mem2`, GATK 4 (`gatk`), `trimmomatic`, `fastqc`, `multiqc`, `mosdepth`
- for ONT input: `minimap2`, `medaka`, `chopper`, `NanoPlot`, `multiqc`, `mosdepth`

If you have only the Debian/Ubuntu `trimmomatic.jar`, pass
`--trimmomatic_jar /usr/share/java/trimmomatic.jar`.

Notes for container runs:
- The Python image has no `openpyxl`, so give `.xlsx` clinical metadata as CSV, or use
  `-profile conda`.
- The FHIR validator downloads its IG packages (about 50 MB) into the task directory on each run.
  Without containers it caches them in `~/.fhir`.

## Quick start

```bash
git clone https://github.com/oucru-id/tb-to-fhir-full.git
cd tb-to-fhir-full

# check the installation with the bundled test data
nextflow run . -profile test_stub -stub-run         # wiring only, no tools needed
nextflow run . -profile test_vcf,docker             # end to end on synthetic VCFs
nextflow run . -profile test,docker                 # Illumina path on a 10 MB public subset (ERR1873531, MDR-TB)
```

## Input

### Samplesheet (recommended)

```csv
sample,platform,fastq_1,fastq_2,vcf
TB0001,illumina,/data/TB0001_R1.fastq.gz,/data/TB0001_R2.fastq.gz,
TB0002,ont,/data/TB0002.fastq.gz,,
TB0003,vcf,,,/data/TB0003.vcf.gz
```

```bash
nextflow run oucru-id/tb-to-fhir-full -r v2.0.0 --input samplesheet.csv --outdir results
```

Sample ids may use letters, digits, `.`, `_` and `-` (up to 40 characters). Use laboratory
accessions.

### Folders

```bash
nextflow run . --reads_dir /data/NGS --vcf_dir /data/VCF
```

The folder files must be named:
- `<id>_1_illumina.fastq.gz` and `<id>_2_illumina.fastq.gz`
- `<id>_ont.fastq.gz`
- `<id>.vcf` or `<id>.vcf.gz`

ONT input also needs `--medaka_model`, for example, `r1041_e82_400bps_sup_variant_v5.0.0`.

> **VCF input cannot produce susceptible calls.** A variant-only VCF records which variants were
> *called*, never which regions were *callable*. Coverage therefore cannot be established, and
> every non-resistant drug is reported **Indeterminate**. Resistance detection is unaffected. Use
> FASTQ input if you need susceptibility results.

## FHIR output

The output has three layers:

| Step | Output | Contains |
|---|---|---|
| genomic (always) | `fhir/genomic/<id>.genomic.fhir.json` | DiagnosticReport, variant / susceptibility-panel / region-studied / lineage Observations, Device, Provenance. **No patient data**: the subject is the reference `Patient/<id>-patient`. |
| clinical (with `--patient_metadata`) | `fhir/clinical/` | Patient, Specimen, Encounter, ServiceRequest per sample; Organization, Practitioner, PractitionerRole once. **Personal data.** |
| link (with `--patient_metadata`) | `fhir/bundle/<id>.bundle.fhir.json` | One FHIR `transaction` per sample, with every reference resolved inside the bundle. **Personal data.** |

```bash
nextflow run . --input samplesheet.csv \
  --patient_metadata site/patients.csv \
  --organization_metadata site/organization.csv \
  --practitioner_metadata site/practitioner.csv \
  --performer_org_id 100000001
```

The CSV columns are described in [docs/clinical_metadata.md](docs/clinical_metadata.md). Keep
these files, and the `fhir/clinical` and `fhir/bundle` outputs, out of version control.

A sample with no row in the patient table gets a genomic bundle only. It is listed as
`missing_clinical` in `fhir/bundle/link_summary.tsv`

### Profiles

- `--fhir_profile generic` (default) emits only HL7 Genomics Reporting IG 3.0.0 content.
- `--fhir_profile satusehat` adds the Kemkes StructureDefinitions and `genomics` tag, NIK and
  `sys-ids.kemkes.go.id` identifiers, the administrative-code and citizenship extensions, and the
  Kemkes WHO-grade and conclusion codings.

### Validation

Every bundle is checked with the HL7 FHIR validator against `hl7.fhir.uv.genomics-reporting#3.0.0`
and the project CodeSystems in `assets/fhir/`.
- Results: `fhir/validation/validation_summary.tsv`.

### Upload

Upload runs only when `--fhir_server_url` is given. It sends only linked bundles that passed
validation. Credentials are Nextflow secrets, never parameters:

```bash
nextflow secrets set TBFHIR_CLIENT_ID     '<client id>'
nextflow secrets set TBFHIR_CLIENT_SECRET '<client secret>'
nextflow run . --input samplesheet.csv --patient_metadata site/patients.csv \
  --fhir_profile satusehat \
  --fhir_server_url https://<server>/fhir --fhir_auth_url https://<sso>/realms/<realm>
```

To also send an `X-Api-Key` header, set the `TBFHIR_API_KEY` secret and add
`--fhir_upload_api_key`. For SATUSEHAT, test against the sandbox first.
`--bundle_request_method POST` sends server-assigned-id transactions (with `ifNoneExist` on NIK).
`--external_refs` maps samples to Patients that already exist on the server.

## Configuration

All parameters are described in [nextflow_schema.json](nextflow_schema.json) and listed by
`nextflow run . --help`. The most used:

| Parameter | Default | Purpose |
|---|---|---|
| `--input` | | samplesheet |
| `--outdir` | `results` | output folder |
| `--fhir_profile` | `generic` | `generic` or `satusehat` |
| `--patient_metadata` | | enables the clinical and link steps |
| `--performer_org_id` | | performing laboratory, as the genomic report's performer |
| `--coverage_min_depth` / `--coverage_min_breadth` | `10` / `0.95` | assessability threshold per locus |
| `--filter_min_depth` / `--filter_min_quality` | `5` / `20` | variant call filters |
| `--max_cpus` / `--max_memory` / `--max_time` | `8` / `16 GB` / `24h` | cap on any single task |

### WHO annotation table

The table must have **one row per allele**, with positionally aligned `DRUG` /
`WHO_CLASSIFICATION` / `VARIANT_ID` arrays. A table with one row per (variant, drug) pair causes
`bcftools annotate` to collapse the grade, which silently drops resistance calls. The pipeline
refuses to run on such a table.

```bash
tools/annotation/rebuild_annotation_table.sh      # builds, compresses, indexes and verifies
nextflow run . --annotation_table data/<prefix>_annotation_table.tsv.gz \
               --annotation_header data/<prefix>_annotations_header.txt
```

## Output structure

```
results/
├── pipeline_info/          execution report, timeline, trace, DAG, software_versions.yml
├── qc/fastqc/              FastQC reports (Illumina)
├── qc/nanoplot/            <id>_NanoStats.txt, <id>_NanoPlot-report.html (ONT)
├── qc/multiqc/             multiqc_report.html
├── qc/coverage/            who_targets.bed, <id>.coverage.json, mosdepth beds
├── variants/               <id>.filtered_variants.vcf.gz, <id>.annotated_variants.vcf.gz
├── lineage/                <id>.lineage.json
├── reports/                <id>.summary_report.txt
└── fhir/
    ├── genomic/            <id>.genomic.fhir.json           (PHI-free)
    ├── clinical/           clinical resources, clinical_summary.tsv   (personal data)
    ├── bundle/             <id>.bundle.fhir.json, link_summary.tsv    (personal data)
    ├── validation/         validation_summary.tsv, validation_<n>.json (validator output per batch)
    └── upload/             <id>.bundle.upload.json
```

## Repository layout

```
main.nf  workflows/tbfhir.nf  subworkflows/local/  modules/local/   Nextflow code
bin/                       Python steps (genomic, clinical, link, coverage, lineage, upload)
conf/                      resources, output layout, test profiles
assets/                    samplesheet schema, FHIR CodeSystems, metadata templates
data/                      reference genome, WHO annotation table, lineage barcode, masks
tools/annotation/          rebuild the WHO annotation table
tools/metadata_ui/         Streamlit editor for the clinical metadata table
tools/satusehat/           SATUSEHAT token and patient-lookup helpers
tests/                     pytest suite, synthetic fixtures, v1 equivalence goldens
```

## Support

[GitHub Issues](https://github.com/oucru-id/tb-to-fhir-full/issues)
