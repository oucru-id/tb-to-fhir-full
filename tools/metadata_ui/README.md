# Clinical Metadata UI

Interactive form for the per-sample patient metadata the TB pipeline consumes,
replacing hand-editing of `data/patient_clinical_metadata.csv`.

## Run

```bash
~/medaka-env/bin/streamlit run tools/metadata_ui/app.py
```

Nothing to install — `medaka-env` already has Streamlit, and `sqlite3` is stdlib.
Note the pipeline's own `venv/` does **not** have pandas or Streamlit; use `medaka-env`.

## How it works

- Edits are held in **SQLite** at `data/clinical_metadata.db`.
- On first run the database is seeded from the existing CSV.
- **Export to pipeline CSV** writes `data/patient_clinical_metadata.csv`, which is
  what `scripts/clinical_metadata_parser.py` reads. Nothing is written to the CSV
  until you press that button.

Every column is stored as `TEXT`. This is deliberate: `nik` and the administrative
codes are all-digit strings whose leading zeros must survive
(`0000012345678901`, not `12345678901`).

## Sample discovery

Sample IDs are discovered the same way `main.nf` builds its channels:

| Input | Glob | Sample ID |
|---|---|---|
| Illumina | `data/NGS/*_{1,2}_illumina.fastq.gz` | strip the pair suffix |
| Nanopore | `data/NGS/*_ont.fastq.gz` | strip `_ont.fastq.gz` |
| VCF | `data/VCF/*.vcf`, `*.vcf.gz` | strip the extension |

Nanopore IDs are stored **without** the `_ont` suffix; `find_matching_sample()`
normalises it away when the pipeline looks the row up.

Samples are shown in three states:

- **Needs metadata** — an input exists but no row. The pipeline silently skips these:
  no Patient, Specimen or ServiceRequest is created, while the DiagnosticReport still
  references a Patient that will not exist.
- **Complete** — has a row.
- **Orphaned** — a row exists but the input file is gone.

## CSV contract

Written by `export.py`, and the format must not drift:

- Exactly the 29 columns in `schema.COLUMNS`, in that order. `updated_at` is internal
  and is not exported.
- Plain UTF-8, no BOM.
- `csv.DictWriter`, so addresses containing commas are quoted correctly.
- Empty values are written as `""`, never `nan` or `None`.
- Written to a temp file then `os.replace()`d, so a crash cannot truncate the file the
  pipeline is about to read.
- The previous file is copied to `data/patient_clinical_metadata.csv.bak` first.

## Out of scope

Organization and practitioner metadata are single-record config files — both loaders
only ever read row 0 — so they are still edited by hand in
`data/organization_metadata.csv` and `data/practitioner_metadata.csv`.

`age`, `diagnosis_date`, `treatment_history` and `contact_number` are collected here
but are not currently mapped into FHIR by anything in `scripts/`.

## Files

| File | Purpose |
|---|---|
| `schema.py` | The 29-column contract, column widgets, new-row defaults |
| `store.py` | SQLite init / load / upsert / seed |
| `discovery.py` | Scan `data/NGS` and `data/VCF` for sample IDs |
| `export.py` | Atomic CSV write with backup, and a diff against the file on disk |
| `app.py` | The Streamlit page |
