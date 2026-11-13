# Changelog

All notable changes to TBtoFHIR are documented here.
Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project uses semantic versioning.

---

## [2.0.0] — Unreleased

Restructures the pipeline for public release:
- The FHIR output is split into a **PHI-free genomic layer** and an **optional clinical layer**.
- Plain HL7 Genomics Reporting is the default profile; SATUSEHAT becomes opt-in.
- The Nextflow code is reorganised into modules with stubs, a parameter schema and test profiles.

With `--fhir_profile satusehat` the linked bundles match 1.6 except for the fixes listed below
(checked by `tests/test_v1_equivalence.py` against `tests/compat/allowlist.py`).

### Fixed

- **WHO grades came from the initial instead of the final grading.** The annotation table was
  built from `Initial confidence grading ALL dataset`, which leaves out the WHO additional grading
  criteria (LoF rule, INH-ETO and BDQ-CFZ cross-resistance, WHO-endorsed gDST assays, selection
  evidence). 925 variant-drug pairs that WHO grades *Assoc w R* or *Assoc w R - Interim* were
  reported as *Uncertain significance*, so the drug was called susceptible. Examples: fabG1
  promoter `inhA_c.-779G>T` and `inhA_c.-770T>A/C/G` (INH, ETH), 757 named frameshift and
  stop-gained variants upgraded under the LoF rule (e.g. `katG_p.Gly124fs`, `ethA_p.Lys200fs`,
  `ddn_p.Trp88*`), `atpE_p.Glu61Asp` (BDQ) and `Rv0678_p.Leu117Arg` (CFZ). Another 12,332 pairs
  graded *Not assoc w R - Interim* were shown as *Uncertain significance*. The table is rebuilt
  from `FINAL CONFIDENCE GRADING` as `data/2026-11-12_annotation_table.tsv.gz` and the 2026-09-07
  table is removed. Allele rows, genes and coverage targets are unchanged.
  `prepare_who_annotation.py` now stops if the final grading column is missing instead of falling
  back to another column.
- **Illumina calls were filtered on genotype quality (GQ).** GATK's diploid GQ is low for
  near-fixed variants with a few reference reads, so real resistance variants were dropped
  (e.g. `gyrA_p.Asp94Gly` at 119/128 reads, QUAL 4247, GQ 15). Illumina calls are now filtered on
  `QUAL >= filter_min_quality`; ONT and VCF input are unchanged.
- **Samples could receive another sample's results.** Lineage, coverage and clinical rows were
  also matched on the text before the first `_`, so `TB_001` and `TB_002` were both keyed `TB`.
  One could be gated by the other's coverage (turning *Indeterminate* into *Susceptible*), share a
  lineage, and get no Patient at all. Matching is now exact, and each FHIR step receives exactly
  one lineage and one coverage file per sample. Sample ids without `_` (e.g. ERR accessions) were
  not affected. **Re-run 1.6 results whose sample ids contain `_` and share a prefix.**
- **Dangling references.** v1 always pointed the DiagnosticReport and ServiceRequest at an
  `Encounter/<id>-encounter` that was never created, and wrote a DiagnosticReport for samples
  without a Patient. An Encounter is now created. Samples without clinical data get no
  transaction bundle. The link step fails on any reference that does not resolve.
- **Placeholders are no longer written.** v1 wrote `"Unknown"` names and birth dates (an invalid
  `birthDate`), the run time as specimen collection/receipt time, a 2 mL default quantity, and
  `unknown-org` / `unknown-practitioner` ids. Missing values are now left out and flagged in
  `clinical_summary.tsv`.
- **Metadata errors were silent.** A malformed table yielded zero patients. Missing columns,
  duplicate ids, unreadable files and ids that collide after FHIR-id sanitising are now errors.
- **Repeat masking was skipped silently** when the mask file path did not exist inside the task.
  The mask is now a staged input; `NO_FILE` disables it explicitly with a warning.
- **ONT sample ids kept `.fastq`** (`X_ont.fastq.gz` gave sample `X_ont.fastq`).
- **Command-line parameters on Nextflow 26.04.** 26.04 passes every command-line value as text.
  Every number given on the command line failed schema validation, `--fhir_validation_batch_size`
  could not batch, and `false` counted as true. This affected `--skip_fhir_validation false`,
  `--require_clinical false`, `--publish_clinical false`, `--fail_on_validation_error false`,
  `--legacy_annotation_table false` and `--fhir_upload_api_key false`. Values are now converted to
  their schema type before validation and where the pipeline reads them. This works the same on
  Nextflow 25.04, 25.10 and 26.04.
- **Validation results of all but one batch were lost.** Files are validated in batches of
  `--fhir_validation_batch_size` (50). Every batch published `validation.json` and
  `validation_summary.tsv` to the same path, so a run with more than 50 FHIR files kept only one
  batch's results. The summaries are now merged into one `validation_summary.tsv`, and each
  batch's validator output is kept as `validation_<n>.json`.
- **VCF-input quality filter** named `INFO/GQ` even when only `FORMAT/GQ` was defined; bcftools
  rejects such expressions. It is now built from the tags present.
- **`filter_min_depth` / `filter_min_quality`** were recorded on the Device but not used. The
  filters had 5 and 20 hard-coded. They now come from the parameters, with the same defaults.
- **g.HGVS** was `g.{pos}{ref}>{alt}` for every allele, which is valid only for SNVs. Deletions,
  insertions and delins now use correct notation. The HGVS system URI is `http://varnomen.hgvs.org`.
- **`tbd-codes-cs` codes that are not in the IG** (`coverage-breadth`, `variant-quality`,
  `mapping-quality`) now use a project CodeSystem, published in `assets/fhir/`.
- The 5′/3′ UTR sequence-ontology effects never matched, because of a case mismatch.
- A variant record that failed to convert now fails the step, instead of producing a report with
  the variant missing.
- Patient names no longer appear outside the Patient resource (DiagnosticReport and Specimen
  displays, Specimen note, HTML report).
- Upload credentials are no longer passed on the command line, which made them visible in
  process listings and `.command.sh`.

### Added

- `tbfhir_genomic.py`: the genomic step. It builds the DiagnosticReport, Observations, Device and
  Provenance, refers to the patient only as `Patient/<id>-patient`, and adds a performer only
  with `--performer_org_id`.
- `clinical_to_fhir.py`: the clinical step. It builds Patient, Specimen, **Encounter** and
  ServiceRequest per sample, plus Organization, Practitioner and PractitionerRole.
- `link_fhir_bundle.py`: the link step. It builds one `transaction` per sample with references
  resolved, PUT or POST (`ifNoneExist` on NIK), and `--external_refs` for resources already on the
  server.
- `--fhir_profile generic|satusehat`. In `generic`:
  - DiagnosticReport uses the IG 3.0.0 `genomic-report` profile and LOINC 51969-4 (81247-9 kept as
    a second coding).
  - WHO grades use the project `who-catalogue-grade` CodeSystem.
  - The WHO variant name is recorded in its own component, because a valid p.HGVS needs a protein
    accession.
  - Generic output validates with 0 errors against IG 3.0.0.
- Samplesheet input (`--input`, validated against `assets/schema_input.json`).
- `nextflow_schema.json` with parameter validation (nf-schema), and a `manifest`.
- Reference indexing once per run (`PREPARE_REFERENCE`). v1 rebuilt the bwa-mem2 index and the
  GATK dictionary in every sample.
- FHIR validation that can fail the run: `--fail_on_validation_error`. It runs offline by default
  (`--fhir_tx_server n/a`), and the IG is pinned to 3.0.0 (was `#current`).
- Upload only when `--fhir_server_url` is set, only for linked bundles that passed validation, one
  at a time. Credentials come from Nextflow secrets.
- `resourceLimits` from `--max_cpus` / `--max_memory` / `--max_time`, retry on out-of-memory, and
  GATK `-Xmx` derived from the task memory. v1 asked for 10 GB, which failed on a 9.7 GB host.
- Profiles `test_stub` (wiring) and `test_vcf` (end to end on synthetic VCFs), and a pytest suite
  with synthetic fixtures and v1 equivalence goldens.
- `software_versions.yml` collected from every process.
- `tools/dev/check_no_phi.py` and a pre-commit configuration (gitleaks).
- **Container and conda profiles**: `docker`, `podman`, `singularity`, `apptainer`, `conda`,
  `mamba`. Every process has a pinned BioContainers image and bioconda/conda-forge package, at
  the tool versions that produced the 1.6 results.
- GitHub Actions CI: lint, PHI guard, gitleaks, pytest, stub run, the `test_vcf` and `test`
  profiles under Docker, and the conda profile.
- `test` profile: one public Illumina isolate (ENA), reduced to the reads over the WHO resistance
  targets and lineage barcode sites, so the full Illumina path runs in CI.

### Changed

- **ONT read QC uses NanoPlot 1.48.0** instead of FastQC. FastQC's per-base modules assume
  short reads of similar length. The reports are in `qc/nanoplot/`, and MultiQC includes the
  NanoPlot statistics. Illumina reads still use FastQC.
- Scripts moved from `scripts/` to `bin/` and run from `PATH`, which is required for containers.
  Offline tools moved to `tools/annotation/` and `tools/satusehat/`.
- `get_access_token.py` and `get_patient_by_nik.py` keep tokens in `$TBFHIR_HOME`
  (default `~/.tbfhir`), not in the repository. Token files are written owner-only.
- The clinical metadata parser no longer needs pandas.
- Trimmomatic is called as `trimmomatic` (or `--trimmomatic_jar`), and GATK as `gatk` from
  `PATH`. Paths into the installation are no longer hard-coded.
- The `hpc` profile is now `slurm`, with no queue set.
- bwa-mem2 runs with `-K 100000000`, so alignments no longer depend on the number of CPUs. On the
  six v1 validation samples, v1 (1 CPU) and 2.0 (8 CPUs) gave identical variant sets and
  classifications. Depth, QUAL and MQ differed at a handful of variants.
- `upload_fhir.py` uses only the Python standard library (`requests` is no longer needed).
- FHIR validation runs in two steps: the validator (Java image) and its summary (Python image).
  In containers the validator caches IG packages in the task directory.
- The metadata editor stores its database and CSV in `$TBFHIR_METADATA_DIR`
  (default `~/.tbfhir/metadata`, owner-only), not in the repository's `data/`. It finds samples
  from `--samplesheet`, `--reads-dir` or `--vcf-dir`.

### Removed

- `annotated_to_fhir.py`, `merge_clinical_fhir.py`, `get_versions.py`, `workflows/*.nf` (v1
  layout).
- Parameters `version`, `illumina_suffix`, `nanopore_suffix`, `max_threads`, `api_key`, and the
  defaults pointing at `data/*_metadata.csv`.
- `data/access_token.json` and `data/input_sso.json` placeholders.

### Breaking Changes

- **Output layout**: see the README.
  - `fhir/` is now `fhir/genomic/`.
  - `fhir_merged/` is now `fhir/bundle/`.
  - `runningstat/` is now `pipeline_info/`.
  - `qc/` is split into `qc/fastqc`, `qc/nanoplot`, `qc/multiqc` and `qc/coverage`.
- **`--results_dir` is now `--outdir`**, with default `./results` (no longer inside the
  installation).
- **Inputs must be given**: `--input samplesheet.csv`, or `--reads_dir` / `--vcf_dir`. The
  defaults no longer point into `data/`.
- **Clinical metadata is opt-in**: give `--patient_metadata` (and the organization and
  practitioner tables) to get linked bundles.
- **Default FHIR profile is `generic`**. Use `--fhir_profile satusehat` for 1.6-compatible
  SATUSEHAT output.
- **Genomic bundles contain the DiagnosticReport, Device and Provenance.** These used to be added
  by the merge step. The performer is set only with `--performer_org_id`; v1 hard-coded a
  facility id.
- **`--medaka_model` is required for ONT input.** The silent R9.4.1 default is gone.
- **Nextflow ≥ 25.04 is required** (the nf-schema plugin needs it).

---

## [1.6.0] — 2026-09-07

Fixes two defects that silently produced incorrect drug-resistance calls, and changes how
susceptibility is reported. **Results from this version are not directly comparable with v1.5.0** —
see [Breaking Changes](#breaking-changes).

### Fixed

- **Multi-drug variants no longer lose their resistance call.** The WHO catalogue grades each
  *(variant, drug)* pair independently, but `WHO_CLASSIFICATION` was declared `Number=1` while
  `DRUG` was `Number=.`, so `bcftools annotate` kept every drug but only one arbitrary grade.
  `annotated_to_fhir.py` then added the comma-joined string (e.g. `"bedaquiline,clofazimine"`) to
  the resistant-drug set, so the per-drug lookup failed for **both** drugs and the panel reported
  them susceptible. Affected **2,500 catalogue alleles, 2,431 of them bedaquiline** — the drug that
  separates pre-XDR from XDR and decides BPaL versus a longer individualized regimen.
- **WHO group 2 is no longer reported as group 1.** Grade standardization matched
  `'Assoc w R' in classification`, which is true for `"2) Assoc w R - Interim"`, promoting interim
  evidence to confirmed resistance. Affected **519 variant-drug pairs**. Groups 4 and 5 were also
  inverted (latent — the catalogue column contains no group 4 entries).
- **WHO group 3 (Uncertain significance) is no longer folded into "Not assoc w R"** in the summary
  report. It is now counted and listed as its own category, with an explicit note that group 3 is
  not a susceptible result.
- **Genes without a verified NCBI GeneID no longer emit an invalid coding.** The gene *symbol* was
  being placed in the `code` slot under the NCBI Gene system (e.g. `code: "mmpL5"`), affecting
  roughly half of the observed genes. They now emit `text` only.
- **Resistance-gene detection read only the last significance component**, so a gene graded
  `Assoc w R` for one drug and `Uncertain` for another was omitted from the report conclusion.
- **`48013-7` reference-sequence component had no `value[x]`**, making it non-conformant.
- **Records that failed to build were silently discarded** by a bare `try/except`; they are now
  counted and reported with their genomic position.
- Variant Observations now carry `valueRange.high`, so indels have a span.
- **Coverage targets no longer span the whole genome for origin-adjacent genes.** *dnaA* (Rv0001)
  sits at the origin of the circular chromosome, so its promoter variants (`c.-13`) wrap to
  position ~4,411,520. Taking one min→max span per gene produced a single **3.9 Mb** "dnaA"
  interval, which would have made isoniazid assessability depend on breadth across most of the
  genome. Variant positions are now clustered per gene (`--max-gap`, default 5000), so dnaA
  correctly yields two loci. Total targets: 72 intervals / 142 kb.
- **Coverage targets are clamped to contig bounds.** Padding a gene near the end of the
  chromosome produced an interval past the contig end, and `mosdepth` aborts hard on that
  (`index 4411533 not in 0 .. 4411532`). `make_who_targets.py` now takes `--reference-fai` and
  clamps; `MAKE_WHO_TARGETS` generates the index if absent.

### Added

- **Coverage assessment** (`workflows/coverage.nf`, `scripts/make_who_targets.py`,
  `scripts/summarize_coverage.py`). Derives per-drug resistance loci from the WHO catalogue
  (72 intervals, 65 genes, 15 drugs), runs `mosdepth` over them, and computes per-drug
  assessability. VCF-only input produces an explicit "coverage not assessed" record rather than
  silence.
- **`region-studied` Observations** recording per-locus mean depth and breadth, so a reviewer can
  see which loci were assessable and why a drug was left indeterminate.
- **`Provenance` and `Device` resources**, linking the report and every Observation to the software
  version that produced them.
- **Per-drug coded components** on each variant: medication assessed (LOINC `51963-7`) paired with
  that drug's own clinical significance (`53037-8`). Wires up `get_drug_snomed_mapping()`, which
  existed in the codebase but was never called; adds delamanid, pretomanid, cycloserine,
  terizidone and PAS.
- **Allele frequency** (LOINC `81258-6`, from `INFO/DP4`, `AF` or `AF1`), plus **QUAL** and **MQ**
  components. QUAL is now parsed from VCF column 6, which the parser previously discarded.
- **Lineage confidence, score and barcode-SNP support** as Observation components. The fabricated
  `http://tb-lineage.org` code system is replaced by a properly namespaced one and retained as a
  secondary coding for traceability.
- **`scripts/rebuild_annotation_table.sh`** — rebuilds, compresses, indexes and verifies the WHO
  annotation table. `--compress-only` skips the spreadsheet parse.
- **Coordinate cross-check.** When WHO's own `Genomic_coordinates_*.vcf.gz` is present, the
  extracted coordinates are verified against it and the build aborts on any disagreement, which
  would indicate two different catalogue releases.
- **`data/gene_ncbi_map.tsv`** — optional overrides for the 45 catalogue genes with no verified
  NCBI GeneID.
- Configurable parameters: `annotation_table`, `annotation_header`, `coverage_min_depth` (10),
  `coverage_min_breadth` (0.95), `promoter_padding` (200), `filter_min_depth`,
  `filter_min_quality`.

### Changed

- **Susceptibility is now a three-state result**: `Resistant` (`LA6676-6`), `Susceptible`
  (`LA24225-7`), `Indeterminate` (`LA9663-1`). A drug is only reported susceptible when every one
  of its resistance loci met the coverage threshold. Resistance calling is unaffected — a detected
  mutation is a positive finding that carries its own evidence.
- **The annotation table is now one row per allele** with positionally aligned
  `DRUG` / `WHO_CLASSIFICATION` / `VARIANT_ID` arrays (`Number=.`), rather than one row per
  (variant, drug) pair.
- **New resistance classifications** replacing `Sensitive`: `No resistance detected`,
  `No resistance detected - partial`, `Indeterminate`, and
  `Drug-resistant - rifampicin not assessable`. The Kemkes `TB-SO` ("drug-sensitive TB") coding is
  only emitted when every drug group was genuinely assessable.
- **RR/MDR/pre-XDR/XDR are only asserted when rifampicin itself was assessable.** An indeterminate
  rifampicin now blocks classification instead of falling through to a lower class.
- **HR-TB requires rifampicin to be assessed and susceptible**, not merely "not resistant".
- Drug-susceptibility results that could not be assessed are named explicitly in
  `DiagnosticReport.conclusion`.
- Illumina and Nanopore sub-workflows now emit their BAM for coverage analysis.
- Annotation table and header paths moved from hardcoded values to `nextflow.config` parameters.
- `Device.deviceName` is `TBtoFHIR`; `Device.property` carries the source repository URL.
- Catalogue extraction uses the spreadsheet's own `chromosome` column instead of assuming
  `NC_000962.3`.

### Removed

- `file_sha256()` helper and its `hashlib` import (obsoleted by the simplified `Device`).
- `--annotation_table` and `--lineage_barcode` arguments from `merge_clinical_fhir.py`.

### Breaking Changes

- **The annotation table must be rebuilt.** Run `scripts/rebuild_annotation_table.sh` and update
  `annotation_table` / `annotation_header` in `nextflow.config`. The pipeline refuses to run
  against a legacy one-row-per-pair table and reports how to fix it;
  `--legacy-annotation-table` opts into the old permissive behaviour, which over-calls resistance.
- **VCF-only input no longer yields susceptible calls.** A variant-only VCF records which variants
  were *called*, never which regions were *callable*, so every non-resistant drug is reported
  `Indeterminate`. v1.5.0 reported these as `Susceptible` with no supporting evidence. Use FASTQ
  input for susceptibility results.
- **`mosdepth` is a new dependency** for FASTQ input.
- Consumers parsing `Observation.text.div` for the associated drug should migrate to the coded
  `51963-7` / `53037-8` component pairs.

### Known Limitations

- Deletion and some LoF graded-variants have no genomic coordinates in the WHO catalogue
  (`pncA_deletion` is named explicitly in WHO's instructions), so they cannot be matched by
  coordinate lookup. Some pyrazinamide resistance is invisible by construction.
- Allele frequency is derived from `INFO/DP4`, `AF` or `AF1`; callers reporting allele depth only
  in `FORMAT/AD` yield no frequency, so heteroresistance cannot be assessed for those inputs.
- 45 WHO catalogue genes have no verified NCBI GeneID and emit text-only codings.
- FHIR validation reports failures but does not fail the run.
- `Specimen.collectedDateTime` / `receivedTime` fall back to run timestamps when the metadata CSV
  columns are empty.

---

## [1.5.0] — earlier

Initial release covering Illumina, Nanopore and VCF input; WHO catalogue annotation; lineage
classification; FHIR R4 bundle generation; clinical metadata merge; and FHIR server upload.
