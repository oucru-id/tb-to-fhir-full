
import csv
import os
import shutil
import tempfile

from schema import COLUMNS


def export_csv(rows, csv_path):

    csv_path = os.path.abspath(csv_path)
    target_dir = os.path.dirname(csv_path)
    os.makedirs(target_dir, exist_ok=True)

    backup_path = None
    if os.path.exists(csv_path):
        backup_path = csv_path + ".bak"
        shutil.copy2(csv_path, backup_path)

    ordered = sorted(rows, key=lambda r: str(r.get("sample_id", "")))

    handle = tempfile.NamedTemporaryFile(
        mode="w", newline="", encoding="utf-8", dir=target_dir,
        prefix=".patient_metadata_", suffix=".tmp", delete=False,
    )
    try:
        with handle:
            writer = csv.DictWriter(handle, fieldnames=COLUMNS, extrasaction="ignore")
            writer.writeheader()
            for row in ordered:
                writer.writerow({c: _clean(row.get(c)) for c in COLUMNS})
        os.replace(handle.name, csv_path)
    except Exception:
        if os.path.exists(handle.name):
            os.unlink(handle.name)
        raise

    return {"path": csv_path, "rows": len(ordered), "backup": backup_path}


def _clean(value):
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in ("nan", "nat", "none"):
        return ""
    return text


def diff_against_csv(rows, csv_path):
    if not os.path.exists(csv_path):
        return [r.get("sample_id", "") for r in rows], [], []

    on_disk = {}
    with open(csv_path, newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle):
            clean = {(k or "").strip(): (v or "").strip() for k, v in raw.items()}
            if clean.get("sample_id"):
                on_disk[clean["sample_id"]] = clean

    new, changed = [], []
    for row in rows:
        sample_id = str(row.get("sample_id", "")).strip()
        if not sample_id:
            continue
        if sample_id not in on_disk:
            new.append(sample_id)
        else:
            before = {c: _clean(on_disk[sample_id].get(c)) for c in COLUMNS}
            after = {c: _clean(row.get(c)) for c in COLUMNS}
            if before != after:
                changed.append(sample_id)

    stored_ids = {str(r.get("sample_id", "")).strip() for r in rows}
    removed = sorted(set(on_disk) - stored_ids)
    return new, changed, removed
