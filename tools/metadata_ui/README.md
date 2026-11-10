# Clinical Metadata UI

A form for the per-sample patient table the pipeline reads with `--patient_metadata`. It
replaces hand-editing the CSV.

## Run

```bash
mamba env create -f envs/metadata_ui.yml          # once (streamlit + pandas)
mamba activate tbfhir-metadata-ui
streamlit run tools/metadata_ui/app.py -- --samplesheet /path/to/samplesheet.csv
```

Any Python environment with Streamlit and pandas works, e.g. `~/medaka-env/bin/streamlit`.

| Option | Environment variable | Default |
|---|---|---|
| `--metadata-dir` | `TBFHIR_METADATA_DIR` | `~/.tbfhir/metadata` |
| `--samplesheet` | `TBFHIR_SAMPLESHEET` | none |
| `--reads-dir`, `--vcf-dir` | `TBFHIR_READS_DIR`, `TBFHIR_VCF_DIR` | none (v1.6 file naming) |

Everything the app writes stays in the metadata directory, which is created owner-only
(`0700`):
- `clinical_metadata.db`
- `patient_clinical_metadata.csv`
- `patient_clinical_metadata.csv.bak`

Never point the metadata directory inside the repository.

## How it works

- Edits are held in **SQLite** (`clinical_metadata.db`).
- On first run the database is seeded from an existing `patient_clinical_metadata.csv` in the
  metadata directory.
- **Export to pipeline CSV** writes `patient_clinical_metadata.csv`. Pass that file to the
  pipeline as `--patient_metadata`. Nothing is written to the CSV until you press the button.

Every column is stored as `TEXT`, deliberately: `nik` and the administrative codes are all-digit
strings whose leading zeros must survive (`0000012345678901`, not `12345678901`).

## Sample states

Sample ids come from the samplesheet (`sample` column) or from the folders, and are shown in
three groups:

- **Needs metadata**: an input exists but there is no row. The pipeline gives these a genomic
  report only and lists them as `missing_clinical`.
- **Complete**: has a row.
- **Orphaned**: a row exists but no input has that sample id.

## CSV contract

Written by `export.py`; the format must not drift:

- Exactly the 29 columns in `schema.COLUMNS`, in that order. `updated_at` is internal and not
  exported.
- Plain UTF-8, no BOM.
- `csv.DictWriter`, so addresses containing commas are quoted correctly.
- Empty values are written as `""`, never `nan` or `None`. The pipeline leaves blank fields out
  of the FHIR resources.
- Written to a temp file and then `os.replace()`d, so a crash cannot truncate the file.
- The previous file is copied to `patient_clinical_metadata.csv.bak` first.

The column meanings are documented in [docs/clinical_metadata.md](../../docs/clinical_metadata.md).

## Out of scope

- Organization and practitioner metadata are single-record files edited by hand. Templates are
  in `assets/metadata_templates/`.
- `age`, `diagnosis_date`, `treatment_history` and `contact_number` are collected here but not
  mapped into FHIR.

## Files

| File | Purpose |
|---|---|
| `schema.py` | The 29-column contract, column widgets, new-row defaults |
| `store.py` | SQLite init / load / upsert / seed |
| `discovery.py` | Sample ids from a samplesheet or from folders |
| `export.py` | Atomic CSV write with backup, and a diff against the file on disk |
| `app.py` | The Streamlit page |
