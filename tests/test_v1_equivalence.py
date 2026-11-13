"""P2 gate: --fhir_profile satusehat reproduces v1.6 except for the allowlist.

The v1 goldens (tests/golden/v1/) were produced by tests/compat/run_v1.py
from the baseline commit. TB_001 / TB_002 are excluded: v1 mixed their data
up (see tests/test_ids.py).
"""
import pytest

from chain import read, run_all
from compat.allowlist import expected_v2, without_encounters
from conftest import GOLDEN, SAMPLES
from normalize import dumps, normalize_bundle

BUGGY_IN_V1 = {"TB_001", "TB_002"}
WITH_CLINICAL = {"RR01", "MDR01", "PREXDR01", "XDR01", "HR01", "INDEL01"}


@pytest.fixture(scope="module")
def satusehat_run(tmp_path_factory):
    return run_all(SAMPLES, tmp_path_factory.mktemp("satusehat"), profile="satusehat")


@pytest.mark.parametrize("sample", sorted(WITH_CLINICAL))
def test_linked_bundle_matches_v1_within_allowlist(sample, satusehat_run):
    actual = normalize_bundle(read(satusehat_run / "bundle" / f"{sample}.bundle.fhir.json"))
    v1 = read(GOLDEN / "v1" / f"{sample}.merged.json")
    expected = expected_v2(v1, f"{sample}.annotated_variants.vcf")
    assert dumps(without_encounters(actual)) == dumps(expected)


@pytest.mark.parametrize("sample", sorted(set(SAMPLES) - BUGGY_IN_V1))
def test_genomic_observations_match_v1(sample, satusehat_run):
    """The PHI-free genomic bundle carries v1's Observations unchanged (allowlist aside)."""
    genomic = normalize_bundle(read(satusehat_run / "genomic" / f"{sample}.genomic.fhir.json"))
    v1 = expected_v2(read(GOLDEN / "v1" / f"{sample}.genomic.json"), f"{sample}.annotated_variants.vcf")
    obs = lambda b: [e for e in b["entry"] if e["resource"]["resourceType"] == "Observation"]  # noqa: E731
    assert dumps({"entry": obs(genomic)}) == dumps({"entry": obs(v1)})


def test_samples_without_clinical_rows_are_not_linked(satusehat_run):
    summary = (satusehat_run / "bundle" / "link_summary.tsv").read_text().splitlines()
    status = dict(line.split("\t")[:2] for line in summary[1:])
    assert status["UNC01"] == "missing_clinical"
    assert status["PZA01"] == "missing_clinical"
    assert not (satusehat_run / "bundle" / "UNC01.bundle.fhir.json").exists()
