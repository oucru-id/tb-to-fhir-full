"""SQLite store for patient clinical metadata.

Every column is TEXT on purpose: nik and the administrative codes are all-digit
strings whose leading zeros must survive. Storing them as INTEGER would reproduce
the pandas coercion bug this store exists to avoid.
"""

import csv
import os
import sqlite3
from datetime import datetime, timezone

from schema import COLUMNS

DB_COLUMNS = COLUMNS + ["updated_at"]


def _connect(db_path):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


def init_db(db_path):
    os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
    cols = ",\n    ".join(
        "{} TEXT PRIMARY KEY".format(c) if c == "sample_id" else "{} TEXT".format(c)
        for c in DB_COLUMNS
    )
    with _connect(db_path) as conn:
        conn.execute(
            "CREATE TABLE IF NOT EXISTS patient_metadata (\n    {}\n)".format(cols)
        )
    return db_path


def load_all(db_path):
    """Return every row as a plain dict, sorted by sample_id."""
    with _connect(db_path) as conn:
        rows = conn.execute(
            "SELECT * FROM patient_metadata ORDER BY sample_id"
        ).fetchall()
    return [{k: (r[k] if r[k] is not None else "") for k in r.keys()} for r in rows]


def upsert_rows(db_path, rows):
    """Insert or replace rows keyed on sample_id. Returns the number written."""
    stamped = datetime.now(timezone.utc).isoformat()
    payload = []
    for row in rows:
        sample_id = str(row.get("sample_id", "")).strip()
        if not sample_id:
            continue
        record = {c: str(row.get(c, "") or "").strip() for c in COLUMNS}
        record["sample_id"] = sample_id
        record["updated_at"] = stamped
        payload.append([record[c] for c in DB_COLUMNS])

    if not payload:
        return 0

    placeholders = ", ".join("?" for _ in DB_COLUMNS)
    column_list = ", ".join(DB_COLUMNS)
    sql = "INSERT OR REPLACE INTO patient_metadata ({}) VALUES ({})".format(
        column_list, placeholders
    )
    with _connect(db_path) as conn:
        conn.executemany(sql, payload)
    return len(payload)


def delete_rows(db_path, sample_ids):
    if not sample_ids:
        return 0
    with _connect(db_path) as conn:
        conn.executemany(
            "DELETE FROM patient_metadata WHERE sample_id = ?",
            [(s,) for s in sample_ids],
        )
    return len(sample_ids)


def seed_from_csv(db_path, csv_path):
    """One-time import of an existing CSV.

    Uses csv.DictReader with utf-8-sig rather than pandas: pd.read_csv without
    dtype=str infers the all-digit nik column as int64 and destroys its leading
    zeros. utf-8-sig also transparently handles the BOM the shipped file carries.
    """
    if not os.path.exists(csv_path):
        return 0

    rows = []
    with open(csv_path, newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle):
            clean = {(k or "").strip(): (v or "").strip() for k, v in raw.items()}
            if clean.get("sample_id"):
                rows.append(clean)

    return upsert_rows(db_path, rows)


def is_empty(db_path):
    with _connect(db_path) as conn:
        count_row = conn.execute("SELECT COUNT(*) FROM patient_metadata").fetchone()
    return count_row[0] == 0
