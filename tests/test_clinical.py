"""Clinical metadata parsing and the clinical resources."""
import importlib.util
import json
import re

import pytest

from chain import read, run_clinical
from clinical_metadata_parser import (PATIENT_COLUMNS, MetadataError, get_clinical_value,
                                      load_clinical_metadata, load_organization_metadata)
from conftest import FIXTURES, REPO

PATIENT_CSV = FIXTURES / "metadata" / "patient_clinical_metadata.csv"
ALL = ["RR01", "MDR01", "PREXDR01", "XDR01", "HR01", "INDEL01", "UNC01", "PZA01", "TB_001", "TB_002"]


def test_columns_match_metadata_ui_schema():
    spec = importlib.util.spec_from_file_location("ui_schema", REPO / "tools" / "metadata_ui" / "schema.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.COLUMNS == PATIENT_COLUMNS


def test_parser_keeps_text_and_returns_none_for_blanks():
    rows = load_clinical_metadata(str(PATIENT_CSV))
    assert rows["RR01"]["nik"] == "9901010000000001"
    assert get_clinical_value(rows["INDEL01"], "birth_date") is None
    assert get_clinical_value(rows["INDEL01"], "birth_date", "fallback") == "fallback"


@pytest.mark.parametrize("content,message", [
    ("family_name\nX\n", "no sample_id column"),
    ("sample_id,family_name\nA,x\nA,y\n", "duplicate sample_id"),
    ("sample_id,family_name\nTB_1,x\nTB.1,y\n", "collide"),
    ("sample_id,family_name\n,x\n", "has no sample_id"),
])
def test_parser_raises_on_bad_input(tmp_path, content, message):
    path = tmp_path / "bad.csv"
    path.write_text(content)
    with pytest.raises(MetadataError, match=message):
        load_clinical_metadata(str(path))


def test_parser_raises_on_missing_file(tmp_path):
    with pytest.raises(MetadataError, match="not found"):
        load_clinical_metadata(str(tmp_path / "nope.csv"))


def test_header_only_organization_file_is_an_error():
    with pytest.raises(MetadataError, match="no data row"):
        load_organization_metadata(str(REPO / "assets" / "metadata_templates" / "organization_metadata.csv"))


@pytest.fixture(scope="module", params=["generic", "satusehat"])
def clinical_run(request, tmp_path_factory):
    out = tmp_path_factory.mktemp(f"clinical-{request.param}")
    run_clinical(ALL, out, profile=request.param)
    return out


def resources(path):
    return {e["resource"]["resourceType"]: e["resource"] for e in read(path)["entry"]}


def test_no_placeholders_anywhere(clinical_run):
    for path in clinical_run.glob("*.json"):
        text = path.read_text(encoding="utf-8")
        assert "Unknown" not in text and "unknown-" not in text, path.name


def test_incomplete_row_omits_missing_elements(clinical_run):
    r = resources(clinical_run / "INDEL01.clinical.json")
    assert "name" not in r["Patient"] and "birthDate" not in r["Patient"]
    assert "receivedTime" not in r["Specimen"]
    assert "collectedDateTime" not in r["Specimen"]["collection"]
    assert "period" not in r["Encounter"]


def test_encounter_is_created_and_referenced(clinical_run):
    r = resources(clinical_run / "RR01.clinical.json")
    enc = r["Encounter"]
    assert enc["id"] == "RR01-encounter" and enc["status"] == "finished"
    assert enc["period"] == {"start": "2026-01-10"}
    assert enc["subject"] == {"reference": "Patient/RR01-patient"}
    assert enc["serviceProvider"] == {"reference": "Organization/99000001"}
    assert r["ServiceRequest"]["encounter"]["reference"] == "Encounter/RR01-encounter"


def test_patient_name_appears_only_in_patient(clinical_run):
    patient_family = re.compile(r"Sintetis\d")    # patients are Sintetis1..9; the practitioner is not
    for path in clinical_run.glob("*.clinical.json"):
        if path.name == "shared.clinical.json":
            continue
        for rtype, resource in resources(path).items():
            if rtype != "Patient":
                assert not patient_family.search(json.dumps(resource)), f"{path.name}:{rtype}"


def test_summary_reports_missing_orphans_and_flags(clinical_run):
    lines = (clinical_run / "clinical_summary.tsv").read_text().splitlines()
    rows = {ln.split("\t")[0]: ln.split("\t") for ln in lines[1:]}
    assert rows["UNC01"][2] == "missing_clinical"
    assert rows["PZA01"][2] == "missing_clinical"
    assert rows["ORPHAN01"][2] == "orphan_row"
    assert "birth_date_missing" in rows["INDEL01"][3]
    assert "Sintetis" not in "\n".join(lines) and "990101" not in "\n".join(lines)


def test_without_practitioner_or_organization(tmp_path):
    run_clinical(["RR01"], tmp_path, profile="generic", organization=False, practitioner=False)
    r = resources(tmp_path / "RR01.clinical.json")
    assert "serviceProvider" not in r["Encounter"]
    assert "requester" not in r["ServiceRequest"] and "collector" not in r["Specimen"]["collection"]
    assert read(tmp_path / "shared.clinical.json")["entry"] == []


def test_encounter_id_column_references_existing_encounter(tmp_path):
    csv_path = tmp_path / "patients.csv"
    csv_path.write_text("sample_id,family_name,encounter_id\nRR01,Sintetis1,ENC-SERVER-1\n")
    ids = tmp_path / "ids.txt"
    ids.write_text("RR01\n")
    from conftest import run_script
    run_script("clinical_to_fhir.py", "--patient_metadata", csv_path, "--sample_ids", ids,
               "--outdir", tmp_path, cwd=tmp_path)
    assert "Encounter" not in resources(tmp_path / "RR01.clinical.json")
    assert read(tmp_path / "external_refs.json") == {"RR01": {"Encounter": "ENC-SERVER-1"}}
