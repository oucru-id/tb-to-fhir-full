# Changelog

All notable changes to TBtoFHIR are documented here.
Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project uses semantic versioning.

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
