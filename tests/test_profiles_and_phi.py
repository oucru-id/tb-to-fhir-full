"""Profile switch and the PHI-free guarantee of the genomic bundles."""
import base64
import csv
import json

import pytest

from chain import read, run_all
from conftest import FIXTURES, SAMPLES

PHI_COLUMNS = ["family_name", "given_name", "birth_date", "nik", "contact_number", "address",
               "latitude", "longitude", "village_code"]


@pytest.fixture(scope="module")
def generic_run(tmp_path_factory):
    return run_all(SAMPLES, tmp_path_factory.mktemp("generic"), profile="generic")


@pytest.fixture(scope="module")
def satusehat_run(tmp_path_factory):
    return run_all(SAMPLES, tmp_path_factory.mktemp("satusehat"), profile="satusehat")


def all_text(path):
    """File text plus decoded presentedForm payloads."""
    text = path.read_text(encoding="utf-8")
    bundle = json.loads(text)
    for entry in bundle.get("entry", []):
        for form in entry["resource"].get("presentedForm", []):
            text += base64.b64decode(form["data"]).decode("utf-8")
    return text


def synthetic_phi_values():
    with open(FIXTURES / "metadata" / "patient_clinical_metadata.csv", encoding="utf-8") as handle:
        lines = [ln for ln in handle if not ln.startswith("#")]
    values = set()
    for row in csv.DictReader(lines):
        values.update(row[c] for c in PHI_COLUMNS if row.get(c))
    with open(FIXTURES / "metadata" / "practitioner_metadata.csv", encoding="utf-8") as handle:
        lines = [ln for ln in handle if not ln.startswith("#")]
    for row in csv.DictReader(lines):
        values.update(row[c] for c in ("nik", "name", "phone", "str_kki_number") if row.get(c))
    return {v for v in values if len(v) > 3}


@pytest.mark.parametrize("profile", ["generic", "satusehat"])
def test_genomic_bundles_carry_no_phi(profile, generic_run, satusehat_run):
    root = generic_run if profile == "generic" else satusehat_run
    phi = synthetic_phi_values()
    assert phi, "fixture PHI list is empty"
    for path in (root / "genomic").glob("*.genomic.fhir.json"):
        text = all_text(path)
        leaked = sorted(v for v in phi if v in text)
        assert not leaked, f"{path.name} contains clinical values: {leaked}"
        types = {e["resource"]["resourceType"] for e in json.loads(path.read_text())["entry"]}
        assert types <= {"DiagnosticReport", "Observation", "Device", "Provenance"}


def test_generic_output_has_no_kemkes_content(generic_run):
    for folder in ("genomic", "clinical", "bundle"):
        for path in (generic_run / folder).glob("*.json"):
            text = path.read_text(encoding="utf-8").lower()
            assert "kemkes" not in text, f"{folder}/{path.name} mentions kemkes"
            assert "sys-ids" not in text


def test_satusehat_output_keeps_kemkes_content(satusehat_run):
    bundle = read(satusehat_run / "bundle" / "RR01.bundle.fhir.json")
    patient = next(e["resource"] for e in bundle["entry"] if e["resource"]["resourceType"] == "Patient")
    assert patient["meta"]["profile"] == ["https://fhir.kemkes.go.id/r4/StructureDefinition/Patient"]
    assert any(i["system"] == "https://fhir.kemkes.go.id/id/nik" for i in patient["identifier"])
    assert bundle["meta"]["profile"] == ["https://fhir.kemkes.go.id/r4/StructureDefinition/Bundle"]


def report_of(bundle):
    return next(e["resource"] for e in bundle["entry"] if e["resource"]["resourceType"] == "DiagnosticReport")


def test_generic_report_uses_ig3_profile_and_code(generic_run):
    report = report_of(read(generic_run / "genomic" / "RR01.genomic.fhir.json"))
    assert report["meta"]["profile"] == [
        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomic-report"]
    assert report["code"]["coding"][0]["code"] == "51969-4"
    assert "tag" not in report.get("meta", {})


def test_satusehat_report_keeps_v1_profile_and_code(satusehat_run):
    report = report_of(read(satusehat_run / "genomic" / "RR01.genomic.fhir.json"))
    assert report["meta"]["profile"] == [
        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomics-report"]
    assert report["code"]["coding"][0]["code"] == "81247-9"


def test_generic_who_grades_use_project_codesystem(generic_run):
    bundle = read(generic_run / "genomic" / "XDR01.genomic.fhir.json")
    systems = set()
    for entry in bundle["entry"]:
        for comp in entry["resource"].get("component", []):
            if comp["code"]["coding"][0]["code"] == "53037-8":
                systems.update(c["system"] for c in comp["valueCodeableConcept"].get("coding", []))
    assert systems == {"http://terminology.spheres.id/CodeSystem/who-catalogue-grade"}


def test_genomic_bundle_without_performer_has_no_organization_reference(tmp_path):
    from chain import run_genomic
    path = run_genomic("RR01", tmp_path, profile="generic", performer=False)
    assert "Organization/" not in path.read_text()
