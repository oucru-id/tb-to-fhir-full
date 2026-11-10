"""Column contract for data/patient_clinical_metadata.csv.

The order of COLUMNS is the CSV contract consumed by
scripts/clinical_metadata_parser.py -- it must not change.
"""

# Exactly the 29 columns the pipeline reads, in order.
COLUMNS = [
    "sample_id",
    "family_name",
    "given_name",
    "gender",
    "birth_date",
    "age",
    "nik",
    "diagnosis_date",
    "treatment_history",
    "contact_number",
    "address",
    "city",
    "state",
    "province_code",
    "city_code",
    "district_code",
    "village_code",
    "latitude",
    "longitude",
    "citizenship_status",
    "specimen_type_code",
    "specimen_type_display",
    "specimen_collection_method_code",
    "specimen_collection_method_display",
    "specimen_collection_method_text",
    "specimen_quantity_value",
    "specimen_quantity_unit",
    "specimen_collected_date",
    "specimen_received_date",
]

# Columns rendered as date pickers. Stored and exported as ISO YYYY-MM-DD strings.
DATE_COLUMNS = [
    "birth_date",
    "diagnosis_date",
    "specimen_collected_date",
    "specimen_received_date",
]

# All-digit identifiers whose leading zeros must survive. Never render these as
# numbers -- a NumberColumn would turn 0000012345678901 into 12345678901.
DIGIT_STRING_COLUMNS = [
    "nik",
    "province_code",
    "city_code",
    "district_code",
    "village_code",
    "contact_number",
    "age",
]

# Values create_patient_resource() knows how to normalise.
GENDER_OPTIONS = ["", "laki-laki", "perempuan"]
CITIZENSHIP_OPTIONS = ["", "WNI", "WNA"]

# Specimen defaults matching the current records and the fallbacks hardcoded in
# scripts/merge_clinical_fhir.py:create_specimen_resource().
SPECIMEN_TYPE_OPTIONS = {
    "119334006": "Sputum specimen",
    "258607008": "Bronchoalveolar lavage fluid specimen",
    "119334006-gastric": "Gastric aspirate specimen",
}
COLLECTION_METHOD_OPTIONS = {
    "KOFFP": "Plate Cough",
    "EXPEC": "Expectorated",
    "BRUSH": "Brush",
}
QUANTITY_UNIT_OPTIONS = ["mL", "L"]

# Applied to a freshly added row so an operator only fills what differs.
NEW_ROW_DEFAULTS = {
    "gender": "",
    "citizenship_status": "WNI",
    "treatment_history": "New case",
    "specimen_type_code": "119334006",
    "specimen_type_display": "Sputum specimen",
    "specimen_collection_method_code": "KOFFP",
    "specimen_collection_method_display": "Plate Cough",
    "specimen_collection_method_text": "Sputum collection",
    "specimen_quantity_value": "2",
    "specimen_quantity_unit": "mL",
}

# Columns copied from the most recently saved row when adding a new one, on the
# assumption that a site processes samples from one catchment area.
INHERITED_COLUMNS = [
    "city",
    "state",
    "province_code",
    "city_code",
    "district_code",
    "village_code",
]

# Collected but never consumed by scripts/ -- surfaced so the operator knows.
UNUSED_BY_PIPELINE = ["age", "diagnosis_date", "treatment_history", "contact_number"]


def blank_row(sample_id="", inherit_from=None):
    """Build a new row dict with defaults and optional inherited location fields."""
    row = {col: "" for col in COLUMNS}
    row.update(NEW_ROW_DEFAULTS)
    if inherit_from:
        for col in INHERITED_COLUMNS:
            value = inherit_from.get(col)
            if value:
                row[col] = value
    row["sample_id"] = sample_id
    return row
