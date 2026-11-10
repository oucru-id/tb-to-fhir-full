
import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import discovery
import export
import schema
import store

def parse_args():
    parser = argparse.ArgumentParser(prog="streamlit run tools/metadata_ui/app.py --")
    parser.add_argument("--metadata-dir", default=os.environ.get(
        "TBFHIR_METADATA_DIR", str(Path.home() / ".tbfhir" / "metadata")))
    parser.add_argument("--samplesheet", default=os.environ.get("TBFHIR_SAMPLESHEET"))
    parser.add_argument("--reads-dir", default=os.environ.get("TBFHIR_READS_DIR"))
    parser.add_argument("--vcf-dir", default=os.environ.get("TBFHIR_VCF_DIR"))
    args, _unknown = parser.parse_known_args(sys.argv[1:])
    return args


ARGS = parse_args()
METADATA_DIR = Path(ARGS.metadata_dir).expanduser()
METADATA_DIR.mkdir(parents=True, exist_ok=True)
os.chmod(METADATA_DIR, 0o700)
DB_PATH = METADATA_DIR / "clinical_metadata.db"
CSV_PATH = METADATA_DIR / "patient_clinical_metadata.csv"

DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")

st.set_page_config(page_title="TB Clinical Metadata", layout="wide")


def bootstrap():
    store.init_db(str(DB_PATH))
    if store.is_empty(str(DB_PATH)) and CSV_PATH.exists():
        seeded = store.seed_from_csv(str(DB_PATH), str(CSV_PATH))
        if seeded:
            st.toast("Imported {} existing row(s) from the CSV".format(seeded))


def to_frame(rows):
    frame = pd.DataFrame(rows, columns=schema.COLUMNS)
    if frame.empty:
        frame = pd.DataFrame({c: pd.Series(dtype="object") for c in schema.COLUMNS})
    for col in schema.COLUMNS:
        if col in schema.DATE_COLUMNS:
            frame[col] = pd.to_datetime(frame[col], errors="coerce")
        else:
            frame[col] = frame[col].fillna("").astype(str)
    return frame


def from_frame(frame):
    rows = []
    for record in frame.to_dict(orient="records"):
        row = {}
        for col in schema.COLUMNS:
            value = record.get(col)
            if col in schema.DATE_COLUMNS:
                if pd.isna(value):
                    row[col] = ""
                else:
                    row[col] = pd.Timestamp(value).strftime("%Y-%m-%d")
            else:
                row[col] = "" if value is None else str(value).strip()
        rows.append(row)
    return rows


def column_config():
    text = st.column_config.TextColumn
    select = st.column_config.SelectboxColumn

    cfg = {}
    for col in schema.COLUMNS:
        cfg[col] = text(col.replace("_", " ").title())

    cfg["sample_id"] = text(
        "Sample ID", disabled=True,
        help="Matches the input filename; the primary key",
    )
    cfg["gender"] = select(
        "Gender", options=schema.GENDER_OPTIONS,
        help="Normalised to FHIR male/female by create_patient_resource()",
    )
    cfg["citizenship_status"] = select(
        "Citizenship", options=schema.CITIZENSHIP_OPTIONS
    )
    cfg["nik"] = text(
        "NIK", max_chars=16,
        help="16 digits. Kept as text so leading zeros survive.",
    )
    cfg["address"] = text("Address", width="medium")
    cfg["specimen_type_code"] = select(
        "Specimen code", options=[""] + list(schema.SPECIMEN_TYPE_OPTIONS),
    )
    cfg["specimen_collection_method_code"] = select(
        "Collection code", options=[""] + list(schema.COLLECTION_METHOD_OPTIONS),
    )
    cfg["specimen_quantity_unit"] = select(
        "Qty unit", options=[""] + schema.QUANTITY_UNIT_OPTIONS,
    )

    for col in schema.DATE_COLUMNS:
        cfg[col] = st.column_config.DateColumn(
            col.replace("_", " ").title(), format="YYYY-MM-DD",
        )


    unused_help = "Collected here but not currently mapped into FHIR by bin/"
    for col in schema.UNUSED_BY_PIPELINE:
        label = col.replace("_", " ").title()
        if col in schema.DATE_COLUMNS:
            cfg[col] = st.column_config.DateColumn(
                label, format="YYYY-MM-DD", help=unused_help,
            )
        else:
            cfg[col] = text(label, help=unused_help)

    return cfg


def validate(rows):
    warnings = []
    for row in rows:
        sample_id = row.get("sample_id") or "(no sample_id)"

        nik = row.get("nik", "")
        if nik and (not nik.isdigit() or len(nik) != 16):
            warnings.append(
                "{}: NIK is {} char(s), expected 16 digits".format(sample_id, len(nik))
            )

        for col in schema.DATE_COLUMNS:
            value = row.get(col, "")
            if value and not DATE_RE.match(value):
                warnings.append(
                    "{}: {} is not YYYY-MM-DD ({})".format(sample_id, col, value)
                )

        if not row.get("family_name") and not row.get("given_name"):
            warnings.append("{}: no patient name".format(sample_id))

        for col in ("latitude", "longitude"):
            value = row.get(col, "")
            if value:
                try:
                    float(value)
                except ValueError:
                    warnings.append(
                        "{}: {} is not numeric ({})".format(sample_id, col, value)
                    )
    return warnings


bootstrap()

st.title("TB Clinical Metadata")
st.caption(
    "Editing {} via {}. Changes save to the database; use Export to write the CSV "
    "you pass to the pipeline as --patient_metadata.".format(CSV_PATH, DB_PATH.name)
)

rows = store.load_all(str(DB_PATH))
existing_ids = [r["sample_id"] for r in rows]
found = discovery.discover_sample_ids(ARGS.samplesheet, ARGS.reads_dir, ARGS.vcf_dir)
groups = discovery.group_samples(found, existing_ids)

needs = groups["needs_metadata"]
orphaned = groups["orphaned"]

if not (ARGS.samplesheet or ARGS.reads_dir or ARGS.vcf_dir):
    st.info("No samples to match against: start the app with "
            "`-- --samplesheet samplesheet.csv` (or --reads-dir / --vcf-dir).")

col_a, col_b, col_c = st.columns(3)
col_a.metric("Samples found", len(found))
col_b.metric("With metadata", len(groups["complete"]))
col_c.metric("Needs metadata", len(needs))

if needs:
    with st.container(border=True):
        st.warning(
            "**{} sample(s) have no metadata row.** The pipeline gives them a genomic "
            "report only (no Patient, Specimen, Encounter or linked bundle) and lists "
            "them as missing_clinical.".format(len(needs))
        )
        st.write(", ".join("{} ({})".format(s, found[s]) for s in needs))
        if st.button("Add rows for these samples", type="primary"):
            latest = max(rows, key=lambda r: r.get("updated_at", ""), default=None)
            store.upsert_rows(
                str(DB_PATH),
                [schema.blank_row(s, inherit_from=latest) for s in needs],
            )
            st.rerun()

if orphaned:
    st.info(
        "Rows with no matching input file (renamed or removed samples): "
        + ", ".join(orphaned)
    )

st.subheader("Patient metadata")

if not rows:
    st.info("No rows yet. Point the app at your samples, then use the button above.")
    st.stop()

edited = st.data_editor(
    to_frame(rows),
    column_config=column_config(),
    num_rows="fixed",
    width="stretch",
    height=420,
    key="grid",
)

edited_rows = from_frame(edited)
warnings = validate(edited_rows)
new, changed, removed = export.diff_against_csv(edited_rows, str(CSV_PATH))

left, right = st.columns([1, 2])

with left:
    if st.button("Save to database", type="primary", width="stretch"):
        written = store.upsert_rows(str(DB_PATH), edited_rows)
        st.success("Saved {} row(s)".format(written))
        st.rerun()

    summary = "Versus the CSV on disk: {} row(s), {} new, {} changed".format(
        len(edited_rows), len(new), len(changed)
    )
    if removed:
        summary += ", {} only in the CSV".format(len(removed))
    st.caption(summary)

    if st.button("Export to pipeline CSV", width="stretch"):
        store.upsert_rows(str(DB_PATH), edited_rows)
        result = export.export_csv(store.load_all(str(DB_PATH)), str(CSV_PATH))
        st.success(
            "Wrote {} row(s) to {}".format(result["rows"], Path(result["path"]).name)
        )
        if result["backup"]:
            st.caption(
                "Previous file backed up to {}".format(Path(result["backup"]).name)
            )

with right:
    if warnings:
        with st.expander(
            "{} validation warning(s)".format(len(warnings)), expanded=True
        ):
            for message in warnings:
                st.write("- " + message)
        st.caption("Warnings do not block saving or export.")
    else:
        st.success("No validation warnings")

with st.expander("Notes on fields"):
    st.markdown(
        "- **NIK** and the administrative codes are stored as text. Entering them "
        "as numbers elsewhere loses leading zeros (0000012345678901 becomes "
        "12345678901).\n"
        "- **Blank cells are left out** of the FHIR resources (no placeholder, no run "
        "timestamp), so fill specimen_collected_date / specimen_received_date where known.\n"
        "- **age, diagnosis_date, treatment_history, contact_number** are collected "
        "but not currently mapped into FHIR by anything in bin/.\n"
        "- **Organization and practitioner** metadata are single-record files edited by "
        "hand; templates are in assets/metadata_templates/."
    )
