#!/usr/bin/env python3

import csv
import os

from tbfhir_common import lookup_sample, sanitize_id

PATIENT_COLUMNS = [
    "sample_id", "family_name", "given_name", "gender", "birth_date", "age", "nik",
    "diagnosis_date", "treatment_history", "contact_number", "address", "city", "state",
    "province_code", "city_code", "district_code", "village_code", "latitude", "longitude",
    "citizenship_status", "specimen_type_code", "specimen_type_display",
    "specimen_collection_method_code", "specimen_collection_method_display",
    "specimen_collection_method_text", "specimen_quantity_value", "specimen_quantity_unit",
    "specimen_collected_date", "specimen_received_date",
]

PATIENT_OPTIONAL_COLUMNS = ["encounter_id", "location_id", "country"]

ORGANIZATION_COLUMNS = [
    "org_id", "name", "alias", "type_code", "type_display", "type_text", "phone", "email",
    "address_line", "city", "state", "country", "province_code", "city_code", "district_code",
    "latitude", "longitude",
]
PRACTITIONER_COLUMNS = [
    "practitioner_id", "nik", "name", "phone", "gender", "birth_date", "str_kki_number",
    "qualification_period_start", "role_id", "role_code", "role_display", "role_text",
]

_SAMPLE_ID_ALIASES = ("sample_id", "sampleid")


class MetadataError(Exception):
    pass


def _normalise_header(name):
    return str(name or "").strip().lower().replace(" ", "_")


def _clean(value):
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        lines = list(handle)
    while lines and (lines[0].lstrip().startswith("#") or not lines[0].strip()):
        lines.pop(0)
    if not lines:
        return [], []
    reader = csv.DictReader(lines)
    header = [_normalise_header(h) for h in (reader.fieldnames or [])]
    rows = []
    for raw in reader:
        row = {_normalise_header(k): _clean(v) for k, v in raw.items() if k is not None}
        if any(v for v in row.values()):
            rows.append(row)
    return header, rows


def _read_xlsx(path):
    try:
        import openpyxl
    except ImportError as exc:
        raise MetadataError(f"{path}: reading .xlsx needs openpyxl, which the container image does not "
                            f"include; save the table as CSV, or run with -profile conda") from exc
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    sheet = workbook.worksheets[0]
    header = None
    rows = []
    for values in sheet.iter_rows(values_only=True):
        if header is None:
            first = values[0] if values else None
            if first is None or str(first).lstrip().startswith("#"):
                continue
            header = [_normalise_header(v) for v in values]
            continue
        if not any(v not in (None, "") for v in values):
            continue
        row = {}
        for key, value in zip(header, values):
            if isinstance(value, float) and value.is_integer():
                value = int(value)     # numeric cell 2.0 -> "2" (leading zeros were already lost)
            row[key] = _clean(value)
        rows.append(row)
    return header or [], rows


def read_table(path):
    if not path or not os.path.exists(path):
        raise MetadataError(f"metadata file not found: {path}")
    lower = path.lower()
    try:
        if lower.endswith(".xlsx"):
            return _read_xlsx(path)
        if lower.endswith(".csv"):
            return _read_csv(path)
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise MetadataError(f"{path}: cannot read ({exc})") from exc
    raise MetadataError(f"{path}: unsupported format (use .csv or .xlsx)")


def load_clinical_metadata(path):
    header, rows = read_table(path)
    key = next((k for k in _SAMPLE_ID_ALIASES if k in header), None)
    if key is None:
        raise MetadataError(f"{path}: no sample_id column (columns: {', '.join(header)})")
    known = set(PATIENT_COLUMNS) | set(PATIENT_OPTIONAL_COLUMNS) | set(_SAMPLE_ID_ALIASES)
    unknown = [h for h in header if h and h not in known]
    if unknown:
        print(f"WARNING: {path}: ignoring unknown column(s): {', '.join(unknown)}")
    records = {}
    seen_fid = {}
    for lineno, row in enumerate(rows, 2):
        sample_id = row.get(key)
        if not sample_id:
            raise MetadataError(f"{path}: data row {lineno - 1} has no sample_id")
        if sample_id in records:
            raise MetadataError(f"{path}: duplicate sample_id {sample_id!r}")
        fid = sanitize_id(sample_id)
        if fid in seen_fid:
            raise MetadataError(f"{path}: sample ids {seen_fid[fid]!r} and {sample_id!r} "
                                f"collide as FHIR id {fid!r}")
        seen_fid[fid] = sample_id
        row["sample_id"] = sample_id
        records[sample_id] = row
    return records


def _load_single_row(path, id_column, kind):
    header, rows = read_table(path)
    if id_column not in header:
        raise MetadataError(f"{path}: no {id_column} column in {kind} metadata")
    if not rows:
        raise MetadataError(f"{path}: {kind} metadata has a header but no data row")
    if len(rows) > 1:
        print(f"WARNING: {path}: {len(rows)} {kind} rows; only the first is used")
    if not rows[0].get(id_column):
        raise MetadataError(f"{path}: first {kind} row has no {id_column}")
    return rows[0]


def load_organization_metadata(path):
    return _load_single_row(path, "org_id", "organization")


def load_practitioner_metadata(path):
    return _load_single_row(path, "practitioner_id", "practitioner")


def find_matching_sample(sample_id, metadata_dict):
    """Exact match only (after platform-suffix stripping and FHIR-id sanitising)."""
    return lookup_sample(metadata_dict, sample_id)


def get_clinical_value(metadata, field, default=None):
    if not metadata:
        return default
    value = metadata.get(_normalise_header(field))
    return value if value not in (None, "") else default
