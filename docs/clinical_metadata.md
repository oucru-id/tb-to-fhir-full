# Clinical metadata

The clinical step (`--patient_metadata`, optional) turns three tables into FHIR resources. Their
outputs are kept apart from the PHI-free genomic bundles:
- **patient table**: one row per sample, giving Patient, Specimen, Encounter and ServiceRequest.
- **organization table**: the performing laboratory, giving Organization.
- **practitioner table**: the responsible clinician, giving Practitioner and PractitionerRole.

These files hold personal data. Keep them outside the repository, for example:

```bash
nextflow run oucru-id/tb-to-fhir-full -params-file ~/site/tbfhir.yaml
```

Header-only templates are in `assets/metadata_templates/`. The Streamlit editor in
`tools/metadata_ui/` writes the patient table for you.

## File format

- CSV (UTF-8, a BOM is fine) or `.xlsx` (needs `openpyxl`).
- Every value is read as text, so NIKs and administrative codes keep their leading zeros. In
  Excel, format those columns as Text; Excel itself drops leading zeros from numeric cells.
- Lines starting with `#` before the header are comments.
- A blank cell means "not known". The element is **left out** of the FHIR resource; nothing is
  invented. (v1.6 wrote `Unknown`, the run time, or a 2 mL default instead.)
- Problems are errors, not silent empty output:
  - a missing file
  - a missing `sample_id` column
  - a duplicate `sample_id`
  - two sample ids that become the same FHIR id (`TB_1` and `TB.1` both become `TB-1`)

## Patient table

`sample_id` is required. All other columns are optional. Unknown columns are ignored with a
warning.

| Column | Used for | Notes |
|---|---|---|
| `sample_id` | links the row to the genomic sample | exact match (or the id without an `_ont` / `_illumina` suffix); never a prefix match |
| `family_name`, `given_name` | Patient.name | |
| `gender` | Patient.gender | `laki-laki`/`pria`/`male`/`m`, `perempuan`/`wanita`/`female`/`f` |
| `birth_date` | Patient.birthDate | `YYYY`, `YYYY-MM` or `YYYY-MM-DD`; anything else is left out and flagged |
| `age`, `diagnosis_date`, `treatment_history`, `contact_number` | not used yet | kept for the metadata editor |
| `nik` | Patient.identifier (satusehat only) | 16 digits; anything else is left out and flagged |
| `address`, `city`, `state` | Patient.address | |
| `country` | Patient.address.country | optional (2.0); satusehat defaults to `ID` |
| `province_code`, `city_code`, `district_code`, `village_code` | Kemkes administrativeCode extension (satusehat only) | |
| `latitude`, `longitude` | geolocation extension on the address | both must be numbers |
| `citizenship_status` | Kemkes citizenshipStatus extension (satusehat only) | `WNI`, `WNA` |
| `specimen_type_code`, `specimen_type_display` | Specimen.type (SNOMED CT) | e.g. `119334006` Sputum specimen |
| `specimen_collection_method_code`, `_display`, `_text` | Specimen.collection.method (v2-0488) | e.g. `KOFFP` Plate, Cough |
| `specimen_quantity_value`, `specimen_quantity_unit` | Specimen.collection.quantity | e.g. `2`, `mL` |
| `specimen_collected_date` | Specimen.collection.collectedDateTime, Encounter.period.start | date or date-time with a time zone |
| `specimen_received_date` | Specimen.receivedTime | |
| `encounter_id` | reference an Encounter that already exists on the server | optional (2.0); no Encounter is generated |
| `location_id` | Encounter.location (a Location that exists on the server) | optional (2.0) |

## Organization table (first row is used)

| Column | Used for |
|---|---|
| `org_id` (required) | Organization.id; performer of the report and its Observations |
| `name`, `alias` | Organization.name / alias |
| `type_code`, `type_display`, `type_text` | Organization.type (Kemkes organization-type coding in satusehat; text in generic) |
| `phone`, `email` | Organization.telecom |
| `address_line`, `city`, `state`, `country` | Organization.address |
| `province_code`, `city_code`, `district_code` | Kemkes administrativeCode extension (satusehat only) |
| `latitude`, `longitude` | geolocation extension |

## Practitioner table (first row is used)

| Column | Used for |
|---|---|
| `practitioner_id` (required) | Practitioner.id; specimen collector, requester, report performer |
| `name`, `phone`, `gender`, `birth_date` | Practitioner |
| `nik` | Practitioner.identifier (satusehat only) |
| `str_kki_number`, `qualification_period_start` | STR-KKI qualification (satusehat only) |
| `role_id`, `role_code`, `role_display`, `role_text` | PractitionerRole (SNOMED CT role code) |

## What the link step adds

For each sample with a patient row, `link_fhir_bundle.py` writes `fhir/bundle/<id>.bundle.fhir.json`,
a FHIR `transaction`:

- `DiagnosticReport.basedOn` → the ServiceRequest
- `DiagnosticReport.encounter` → the Encounter
- `DiagnosticReport.performer` → the Organization and the Practitioner
- every `Type/id` reference rewritten to the `urn:uuid` of the entry it points to. A reference that
  does not resolve is an error.

`clinical_summary.tsv` and `link_summary.tsv` list each sample's status (`matched`,
`missing_clinical`, `orphan_row`, `linked`) and data-quality flags such as `birth_date_missing`.
They never contain the values themselves.
