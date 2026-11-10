# TB Genomics Mutation Analysis to FHIR Genomics Pipeline (TBtoFHIR)

A Nextflow pipeline for *Mycobacterium tuberculosis* genomic analysis from raw sequencing data or
VCFs. It produces WHO-catalogue drug-resistance calls, lineage and coverage QC, and reports them as
HL7 FHIR R4 resources conforming to the
[Genomics Reporting IG 3.0.0](https://hl7.org/fhir/uv/genomics-reporting/STU3/).

**Version 2.0.0 (unreleased)**. See [CHANGELOG.md](CHANGELOG.md). 2.0 changes the inputs, the
output layout and the FHIR structure. Read the
[breaking changes](CHANGELOG.md#breaking-changes) before upgrading from 1.6.

## Key features

- **Multi-platform input**: Illumina paired-end reads, Oxford Nanopore reads, or VCFs.
- **Per-drug WHO grading**: one variant can be *Assoc w R* for one drug and *Uncertain* for
  another.
- **Coverage-gated susceptibility**: a drug is reported susceptible only when every one of its
  resistance loci was sequenced to the required depth and breadth. Otherwise it is
  **Indeterminate**, never silently susceptible.
- **Lineage classification** from barcode SNPs, with confidence and SNP support.
- **QC**: FastQC/MultiQC and per-locus coverage over the WHO resistance targets.
- **PHI-free genomic FHIR by default**. As in
  [vcf2fhir](https://github.com/elimuinformatics/vcf2fhir), the genomic report refers to the
  patient by id only. Clinical resources are built in a separate, optional step and joined to the
  genomic report in a third step.
- **Two FHIR profiles**:
  - `generic` (default): plain HL7 Genomics Reporting.
  - `satusehat`: adds what the Indonesian Ministry of Health SATUSEHAT platform expects (Kemkes
    profiles, NIK identifiers, Kemkes codings).

## Requirements

You need [Nextflow](https://www.nextflow.io/) ≥ 25.04 and **one** of the following:

| Profile | Provides the tools via |
|---|---|
| `-profile docker` (or `podman`) | pinned [BioContainers](https://biocontainers.pro/) images |
| `-profile singularity` (or `apptainer`) | the same images, e.g. on HPC |
| `-profile conda` (or `mamba`) | pinned bioconda / conda-forge packages |
| none | tools already on `PATH` (see below) |

Every process uses the tool versions that produced the v1.6 results:
- FastQC 0.12.1, Trimmomatic 0.39
- bwa-mem2 2.2.1, minimap2 2.26, samtools 1.19.2
- GATK 4.4.0.0, medaka 2.1.0, chopper 0.12.0
- bcftools 1.19, mosdepth 0.3.6, MultiQC 1.30
- Python 3.12, Java 21

Images are pulled only for the processes that run; for example, Illumina-only runs never
download medaka. The largest images are GATK (~1.5 GB), medaka (2 GB) and chopper (1 GB).

Without a profile you need on `PATH`:
- for all runs: `bcftools`, `samtools`, `tabix`/`bgzip`, Python 3 (standard library only) and
  Java ≥ 17
- for Illumina input: `bwa-mem2`, GATK 4 (`gatk`), `trimmomatic`, `fastqc`, `multiqc`, `mosdepth`
- for ONT input: `minimap2`, `medaka`, `chopper`

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

GitHub Actions (`.github/workflows/ci.yml`) runs on every push and pull request:
- the PHI guard, gitleaks, ruff and `nextflow lint`
- pytest
- the stub run on Nextflow 25.04 and the latest release
- the three test profiles under Docker, and `test_vcf` under conda

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
accessions, not national ID numbers: sample ids appear in the PHI-free genomic output.

### Folders (v1.6 naming)

```bash
nextflow run . --reads_dir /data/NGS --vcf_dir /data/VCF
```

The folder files must be named:
- `<id>_1_illumina.fastq.gz` and `<id>_2_illumina.fastq.gz`
- `<id>_ont.fastq.gz`
- `<id>.vcf` or `<id>.vcf.gz`

ONT input also needs `--medaka_model`, for example `r1041_e82_400bps_sup_variant_v5.0.0`.

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
`missing_clinical` in `fhir/bundle/link_summary.tsv`; `--require_clinical` turns that into an
error.

### Profiles

- `--fhir_profile generic` (default) emits only HL7 Genomics Reporting IG 3.0.0 content.
- `--fhir_profile satusehat` adds the Kemkes StructureDefinitions and `genomics` tag, NIK and
  `sys-ids.kemkes.go.id` identifiers, the administrative-code and citizenship extensions, and the
  Kemkes WHO-grade and conclusion codings. Its output matches v1.6, apart from the documented
  fixes.

### Validation

Every bundle is checked with the HL7 FHIR validator against `hl7.fhir.uv.genomics-reporting#3.0.0`
and the project CodeSystems in `assets/fhir/`.
- Results: `fhir/validation/validation_summary.tsv`.
- Validation runs offline by default (`--fhir_tx_server n/a`), so no codes are sent to a public
  terminology server. Code checks that need a server are reported as not performed.
- A validation error stops the run for `generic` and only warns for `satusehat` (the validator
  cannot resolve the Kemkes profiles). Override with `--fail_on_validation_error true|false`.
- `--skip_fhir_validation` skips the step.

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

## Susceptibility interpretation

Each drug in the susceptibility panel resolves to one of **three** states:

| State | LOINC | Meaning |
|---|---|---|
| **Resistant** | `LA6676-6` | A WHO group 1 or 2 variant was detected |
| **Susceptible** | `LA24225-7` | No such variant **and** every locus for that drug met the coverage threshold |
| **Indeterminate** | `LA9663-1` | No such variant, but coverage of the relevant loci was not confirmed |

Region-studied Observations record each locus's mean depth and breadth, so a reviewer can see why
a drug was left Indeterminate.

## Drug resistance classification

| Classification | Criteria |
|---|---|
| XDR-TB | MDR/RR + fluoroquinolone resistance + Group A drug resistance |
| Pre-XDR-TB | MDR/RR + fluoroquinolone resistance |
| MDR-TB | Resistance to both isoniazid and rifampicin |
| RR-TB | Rifampicin resistance only |
| HR-TB | Isoniazid resistance, rifampicin **assessed and susceptible** |
| Mono-resistant | Single-drug resistance (streptomycin, ethionamide, pyrazinamide, ethambutol, or ciprofloxacin) |
| Drug-resistant | Any other resistance combination |
| Drug-resistant - rifampicin not assessable | Resistance detected, but rifampicin could not be assessed |
| No resistance detected | No resistance variant, and **all** drug groups were assessable |
| No resistance detected - partial | No resistance variant, but some drugs could not be assessed |
| Indeterminate | No resistance variant, and no drug could be assessed |

## Output structure

```
results/
├── pipeline_info/          execution report, timeline, trace, DAG, software_versions.yml
├── qc/fastqc/  qc/multiqc/
├── qc/coverage/            who_targets.bed, <id>.coverage.json, mosdepth beds
├── variants/               <id>.filtered_variants.vcf.gz, <id>.annotated_variants.vcf.gz
├── lineage/                <id>.lineage.json
├── reports/                <id>.summary_report.txt
└── fhir/
    ├── genomic/            <id>.genomic.fhir.json           (PHI-free)
    ├── clinical/           clinical resources, clinical_summary.tsv   (personal data)
    ├── bundle/             <id>.bundle.fhir.json, link_summary.tsv    (personal data)
    ├── validation/         validation.json, validation_summary.tsv
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

Developers can run `pytest`. `tests/compat/` checks that `--fhir_profile satusehat` reproduces
v1.6 output, apart from the documented allowlist in `tests/compat/allowlist.py`.

## Known limitations

- **Deletions and some LoF graded variants have no genomic coordinates** in the WHO catalogue
  (`pncA_deletion` is named explicitly in WHO's instructions). They cannot be matched by
  coordinate lookup, so some pyrazinamide resistance is invisible by construction.
- **Multi-nucleotide variants are dropped by the variant filter** (`bcftools view -v
  snps,indels`). GATK reports SNVs by default, but a submitted VCF with MNVs loses them.
- **Lineage without barcode evidence.** The barcode's lineage-4 alleles equal the H37Rv
  reference, so they never appear in a VCF. When no barcode SNP matches (lineage-4 samples, or
  poor coverage at the barcode sites), the classifier reports the first lineage in the barcode
  file (`lineage4.2.2.1`) with 0 matched SNPs and confidence `low`. Check `matched_snps` and
  `confidence` in `lineage/<id>.lineage.json`, or the lineage Observation's components, before
  relying on the call.
- **Allele frequency depends on the caller.** It is derived from `INFO/DP4`, `AF` or `AF1`.
  Callers that report allele depth only in `FORMAT/AD` yield no frequency, so heteroresistance
  cannot be assessed for those inputs.
- **45 WHO catalogue genes have no verified NCBI GeneID** and are emitted as text only. Add
  verified IDs to `data/gene_ncbi_map.tsv`.
- **SATUSEHAT upload has not been verified against the production platform.** Test against the
  sandbox first.

## Support

[GitHub Issues](https://github.com/oucru-id/tb-to-fhir-full/issues)
