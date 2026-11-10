# TB Genomics Mutation Analysis to FHIR Genomics Pipeline (TBtoFHIR)

A platform-agnostic Nextflow pipeline for *Mycobacterium tuberculosis* genomic analysis from raw sequencing data or pre-annotated VCFs, producing HL7 FHIR R4 genomics bundles (IG v3.0.0). [Full documentation](https://tb-pipeline-docs.readthedocs.io/)

**Current version: 1.6.0** — see [CHANGELOG.md](CHANGELOG.md). v1.6.0 fixes two defects that produced incorrect resistance calls and changes how susceptibility is reported; read [Susceptibility Interpretation](#susceptibility-interpretation) before comparing results with earlier runs.

## Key Features

- **Multi-platform**: Illumina paired-end short reads, Oxford Nanopore (ONT) long reads, and pre-annotated VCF input.
- **Drug Resistance Analysis**: Identifies mutations based on the WHO TB mutation catalogue, graded **per (variant, drug) pair** — one variant can be *Assoc w R* for one drug and *Uncertain* for another.
- **Coverage-gated susceptibility**: A drug is only reported susceptible when the resistance loci were demonstrably sequenced. Otherwise it is **Indeterminate**, never silently susceptible.
- **Lineage Classification**: *M. tuberculosis* lineages using barcode SNPs, with confidence and barcode-SNP support.
- **Quality Control**: Per-sample FastQC reports aggregated to MultiQC, plus per-locus coverage over the WHO resistance targets.
- **FHIR Compliance**: HL7 FHIR R4 bundles with Variant, Drug Susceptibility, Region Studied, Lineage Observations, DiagnosticReport, Device and Provenance resources.
- **Clinical Integration**: Merges genomic results with patient, organization, and practitioner metadata.

## Requirements

| Tool | Required for |
|---|---|
| Nextflow (≥ 23) | All runs |
| `bcftools`, `samtools`, `tabix`/`bgzip` | All runs |
| `mosdepth` | Coverage assessment (FASTQ input only) |
| `bwa-mem2`, `minimap2`, `fastqc` | FASTQ input |
| GATK, `trimmomatic` | Illumina path |
| `medaka`, `chopper` | Nanopore path |
| Python 3 with `pandas` | FHIR generation and clinical merge |

> The pipeline calls `python3` directly. Ensure the interpreter on `PATH` has `pandas` installed,
> or the FHIR steps will fail on import.

## Installation

### Setup

```bash
git clone https://github.com/oucru-id/tb-to-fhir-full.git
cd tb-to-fhir-full

# Install Nextflow
curl -s https://get.nextflow.io | bash

# Verify
nextflow -v
```

### Build the WHO annotation table

The annotation table must have **one row per allele** with positionally aligned
`DRUG` / `WHO_CLASSIFICATION` / `VARIANT_ID` arrays. A table with one row per (variant, drug) pair
causes `bcftools annotate` to collapse the grade and silently drop resistance calls.

```bash
scripts/rebuild_annotation_table.sh          # builds, compresses, indexes and verifies
```

Then point `nextflow.config` at the result:

```groovy
annotation_table  = "$baseDir/data/<prefix>_annotation_table.tsv.gz"
annotation_header = "$baseDir/data/<prefix>_annotations_header.txt"
```

If `data/Genomic_coordinates_*.vcf.gz` (WHO's own coordinate VCF) is present, the build
cross-checks the extracted coordinates against it and aborts on any disagreement.

The pipeline refuses to run against a legacy one-row-per-pair table and tells you to rebuild.

## Directory Structure

```
tb-to-fhir-full
├── main.nf                             # Main workflow
├── nextflow.config                     # Configuration and parameters
├── workflows/
│   ├── illumina.nf                     # Illumina sub-workflow
│   ├── nanopore.nf                     # Nanopore sub-workflow
│   ├── vcf.nf                          # VCF sub-workflow
│   ├── lineage.nf                      # Lineage classification
│   ├── coverage.nf                     # WHO target coverage / assessability
│   ├── fhir.nf                         # FHIR variants generation
│   ├── validate_fhir.nf                # FHIR validation
│   ├── merge_clinical_data.nf          # Clinical metadata merge
│   ├── upload_fhir.nf                  # FHIR server upload
│   ├── report.nf                       # QC and sample report generation
│   └── utils.nf                        # Utility functions
├── scripts/
│   ├── annotated_to_fhir.py            # VCF-to-FHIR converter
│   ├── clinical_metadata_parser.py     # Patient/org/practitioner parser
│   ├── generate_sample_report.py       # Per-sample text report
│   ├── lineage_classifier.py           # SNP-barcode lineage classifier
│   ├── make_who_targets.py             # WHO resistance loci -> BED
│   ├── merge_clinical_fhir.py          # FHIR genomics + clinical data merger
│   ├── prepare_who_annotation.py       # WHO catalogue xlsx -> annotation table
│   ├── rebuild_annotation_table.sh     # Rebuild + index + verify the catalogue
│   ├── summarize_coverage.py           # mosdepth -> per-drug assessability
│   ├── upload_fhir.py                  # FHIR uploader
│   ├── get_access_token.py             # Standalone token fetcher
│   └── get_versions.py                 # Software version collector
├── data/
│   ├── NGS/                            # Input FASTQ files
│   ├── VCF/                            # Input VCF files
│   ├── H37Rv.fasta                     # Reference genome
│   ├── repetitive_regions.bed          # Exclusion regions
│   ├── *_lineage.bed                   # Lineage barcode SNPs
│   ├── *_annotation_table.tsv.gz       # WHO mutation annotation table
│   ├── *_annotations_header.txt        # VCF INFO header for the table
│   ├── Genomic_coordinates_*.vcf.gz    # WHO coordinate VCF (optional cross-check)
│   ├── gene_ncbi_map.tsv               # Optional gene -> NCBI GeneID overrides
│   ├── patient_clinical_metadata.csv   # Patient metadata
│   ├── organization_metadata.csv       # Organization metadata
│   └── practitioner_metadata.csv       # Practitioner metadata
└── tools/
    └── fhir-validator.jar              # HL7 FHIR validator
```

## Input Data

### Illumina Reads
Place paired-end FASTQ files in `data/NGS/`:
```
data/NGS/SAMPLE_1_illumina.fastq.gz
data/NGS/SAMPLE_2_illumina.fastq.gz
```

### Nanopore Reads
Place single-end FASTQ files in `data/NGS/`:
```
data/NGS/SAMPLE_ont.fastq.gz
```

### Pre-annotated VCFs
Place VCF files (`.vcf` or `.vcf.gz`) in `data/VCF/`.

> **VCF input cannot produce susceptible calls.** A variant-only VCF records which variants were
> *called*, never which regions were *callable*, so coverage cannot be established and every
> non-resistant drug is reported **Indeterminate**. Resistance detection is unaffected. Use FASTQ
> input if you need susceptibility results.

## Usage

### Get Access Token (FHIR Upload)

```bash
python scripts/get_access_token.py
```

### Basic Run

```bash
nextflow run main.nf
```

### Run with FHIR Upload

> Get the access token first before running with upload.

```bash
nextflow run main.nf \
  --fhir_server_url "https://<BASE_URL>/fhir"
```

## Configuration

Key parameters in `nextflow.config`:

| Parameter | Default | Purpose |
|---|---|---|
| `annotation_table` | `data/*_annotation_table.tsv.gz` | WHO catalogue, one row per allele |
| `annotation_header` | `data/*_annotations_header.txt` | VCF INFO header for the table |
| `coverage_min_depth` | `10` | Depth threshold for a base to count as covered |
| `coverage_min_breadth` | `0.95` | Fraction of a locus that must reach that depth for the drug to be assessable |
| `promoter_padding` | `200` | Bases added around each gene when building coverage targets |
| `filter_min_depth` | `5` | Variant-calling depth filter |
| `filter_min_quality` | `20` | Variant-calling quality filter |

## Susceptibility Interpretation

Each drug in the susceptibility panel resolves to one of **three** states:

| State | LOINC | Meaning |
|---|---|---|
| **Resistant** | `LA6676-6` | A WHO group 1 or 2 variant was detected |
| **Susceptible** | `LA24225-7` | No such variant **and** every locus for that drug met the coverage threshold |
| **Indeterminate** | `LA9663-1` | No such variant, but coverage of the relevant loci was not confirmed |

Absence of a resistance mutation is only evidence of susceptibility when the relevant loci were
demonstrably sequenced. Resistance is unaffected by coverage — a detected mutation is a positive
finding that carries its own evidence.

`region-studied` Observations record per-locus mean depth and breadth, so a reviewer can see
exactly which loci were assessable and why a drug was left indeterminate.

## Drug Resistance Classification

The `DiagnosticReport` conclusion is derived using the following order:

| Classification | Criteria |
|---|---|
| XDR-TB | MDR/RR + Fluoroquinolone resistance + Group A drug resistance |
| Pre-XDR-TB | MDR/RR + Fluoroquinolone resistance |
| MDR-TB | Resistance to both Isoniazid and Rifampicin |
| RR-TB | Rifampicin resistance only |
| HR-TB | Isoniazid resistance, Rifampicin **assessed and susceptible** |
| Mono-resistant | Single drug resistance (Streptomycin, Ethionamide, Pyrazinamide, Ethambutol, or Ciprofloxacin) |
| Drug-resistant | Any other resistance combination |
| Drug-resistant – rifampicin not assessable | Resistance detected, but Rifampicin could not be assessed, so RR/MDR/pre-XDR/XDR cannot be determined |
| No resistance detected | No resistance variant, and **all** drug groups were assessable |
| No resistance detected – partial | No resistance variant, but some drugs could not be assessed |
| Indeterminate | No resistance variant and no drug could be assessed |

> The classifications that depend on Rifampicin (RR/MDR/pre-XDR/XDR) are only asserted when
> Rifampicin itself was assessable.

## Output Structure

```
results/
├── qc/
│   ├── multiqc_report.html         # Aggregated QC report
│   ├── who_targets.bed             # WHO resistance loci used for coverage
│   ├── *.coverage.json             # Per-drug assessability summary
│   ├── *.regions.bed.gz            # mosdepth per-target mean depth
│   └── *.thresholds.bed.gz         # mosdepth bases at each depth threshold
├── lineage/
│   └── *.lineage.json              # Per-sample lineage results
├── fhir/
│   └── *.fhir.json                 # FHIR genomics bundles
├── fhir_merged/
│   └── *.merged.fhir.json          # FHIR bundles with clinical data
├── fhir_validated/
│   └── *.validation.txt            # FHIR validation results
├── reports/
│   └── *.summary_report.txt        # Per-sample summary reports
├── runningstat/
│   ├── execution.html              # Nextflow execution report
│   ├── timeline.html               # Timeline report
│   └── dag.html                    # Workflow DAG
└── software_versions.yml           # Software version manifest
```

## Known Limitations

- **Deletion and some LoF graded-variants have no genomic coordinates** in the WHO catalogue
  (`pncA_deletion` is named explicitly in WHO's instructions), so they cannot be matched by
  coordinate lookup. Some pyrazinamide resistance is therefore invisible by construction.
- **Allele frequency depends on the caller.** It is derived from `INFO/DP4`, `AF` or `AF1`; callers
  that report allele depth only in `FORMAT/AD` will not yield a frequency, so heteroresistance
  cannot be assessed for those inputs.
- **45 WHO catalogue genes have no verified NCBI GeneID** and are emitted as text-only codings.
  Add verified IDs to `data/gene_ncbi_map.tsv` to code them.
- **FHIR validation failures are not fatal** — `validate_fhir.nf` currently reports but does not
  fail the run.
- `Specimen.collectedDateTime` / `receivedTime` fall back to run timestamps when the metadata CSV
  columns are empty.

## Support

[GitHub Issues](https://github.com/oucru-id/tb-to-fhir-full/issues)
