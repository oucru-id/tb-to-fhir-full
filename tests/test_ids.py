"""Sample identity: no prefix matching (the v1.6 TB_001/TB_002 bug)."""
import pytest

from chain import read, run_all
from conftest import FIXTURES, run_script
from tbfhir_common import fhir_id, lookup_sample


def test_lookup_is_exact():
    records = {"TB_001": "a", "TB_002": "b", "X": "x"}
    assert lookup_sample(records, "TB_001") == "a"
    assert lookup_sample(records, "TB_002") == "b"
    assert lookup_sample(records, "TB") is None
    assert lookup_sample(records, "TB_003") is None


def test_lookup_strips_platform_suffix_and_sanitises():
    assert lookup_sample({"S1": 1}, "S1_ont") == 1
    assert lookup_sample({"S1": 1}, "S1_illumina") == 1
    assert lookup_sample({"S1": 1}, "S1_ont.fastq") == 1
    assert lookup_sample({"TB-001": 1}, "TB_001") == 1


def test_lookup_rejects_ambiguous_sanitised_ids():
    with pytest.raises(ValueError):
        lookup_sample({"TB_001": 1, "TB.001": 2}, "TB-001")


def test_fhir_id_limit():
    assert fhir_id("TB_001", "patient") == "TB-001-patient"
    with pytest.raises(ValueError):
        fhir_id("x" * 70)


def panel_state(bundle, loinc):
    panel = next(e["resource"] for e in bundle["entry"] if e["resource"]["id"].endswith("susceptibility-panel"))
    comp = next(c for c in panel["component"] if c["code"]["coding"][0]["code"] == loinc)
    return comp["valueCodeableConcept"]["coding"][0]["display"]


@pytest.fixture(scope="module")
def tb_run(tmp_path_factory):
    return run_all(["TB_001", "TB_002"], tmp_path_factory.mktemp("tb"), profile="satusehat")


def test_underscore_samples_keep_their_own_results(tb_run):
    g1 = read(tb_run / "genomic" / "TB_001.genomic.fhir.json")
    g2 = read(tb_run / "genomic" / "TB_002.genomic.fhir.json")
    assert panel_state(g1, "89489-9") == "Resistant"        # rifampicin
    assert panel_state(g1, "89488-1") == "Susceptible"      # isoniazid, katG covered
    assert panel_state(g2, "89488-1") == "Indeterminate"    # TB_002's own katG coverage is poor
    lineage2 = next(e["resource"] for e in g2["entry"] if e["resource"]["id"] == "TB-002-lineage")
    assert lineage2["valueCodeableConcept"]["coding"][0]["code"] == "lineage2"


def test_underscore_samples_get_their_own_patient(tb_run):
    b1 = read(tb_run / "bundle" / "TB-001.bundle.fhir.json")
    b2 = read(tb_run / "bundle" / "TB-002.bundle.fhir.json")
    ids1 = {e["resource"]["id"] for e in b1["entry"]}
    ids2 = {e["resource"]["id"] for e in b2["entry"]}
    assert {"TB-001-patient", "TB-001-genomic-report"} <= ids1
    assert {"TB-002-patient", "TB-002-genomic-report"} <= ids2
    assert not any(i.startswith("TB-002") for i in ids1)
    assert not any(i.startswith("TB-001") for i in ids2)


def test_mismatched_coverage_file_is_rejected(tmp_path):
    """A coverage file for another sample must never be applied silently."""
    vcf = FIXTURES / "vcf" / "TB_001.annotated_variants.vcf"
    proc = run_script("tbfhir_genomic.py", "--input", vcf, "--sample_id", "TB_001",
                      "--output", tmp_path / "x.json",
                      "--coverage_json", FIXTURES / "coverage" / "TB_002.coverage.json",
                      cwd=tmp_path, check=False)
    assert proc.returncode != 0
    assert "not 'TB_001'" in proc.stderr
