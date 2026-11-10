
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

DATE_COLUMNS = [
    "birth_date",
    "diagnosis_date",
    "specimen_collected_date",
    "specimen_received_date",
]

DIGIT_STRING_COLUMNS = [
    "nik",
    "province_code",
    "city_code",
    "district_code",
    "village_code",
    "contact_number",
    "age",
]

GENDER_OPTIONS = ["", "laki-laki", "perempuan"]
CITIZENSHIP_OPTIONS = ["", "WNI", "WNA"]

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

INHERITED_COLUMNS = [
    "city",
    "state",
    "province_code",
    "city_code",
    "district_code",
    "village_code",
]

UNUSED_BY_PIPELINE = ["age", "diagnosis_date", "treatment_history", "contact_number"]


def blank_row(sample_id="", inherit_from=None):
    row = {col: "" for col in COLUMNS}
    row.update(NEW_ROW_DEFAULTS)
    if inherit_from:
        for col in INHERITED_COLUMNS:
            value = inherit_from.get(col)
            if value:
                row[col] = value
    row["sample_id"] = sample_id
    return row
