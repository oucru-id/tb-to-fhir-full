import json
import random
from functools import lru_cache

import pytest
from conftest import REPO, run_script
from lineage_classifier import classify_lineage, load_barcode_bed, read_pileup, read_reference

BARCODE = REPO / "data" / "2026-05-11_lineage.bed"
REFERENCE = REPO / "data" / "H37Rv.fasta"


@lru_cache(maxsize=1)
def barcode():
    return load_barcode_bed(str(BARCODE))


@lru_cache(maxsize=1)
def reference():
    return read_reference(str(REFERENCE))


def allele_for(marker, lineages):
    ref = reference()[marker["chrom"]][marker["pos"] - 1]
    if marker["lineage"] in lineages:
        return marker["allele"]
    if marker["allele"] == ref:
        return next(base for base in "ACGT" if base != ref)
    return ref


def isolate(path, depth=30, minor=None, minor_fraction=0.0, extra=None):
    markers, _ = barcode()
    evidence = {}
    for marker in markers:
        major = allele_for(marker, set(path) | ({marker["lineage"]} if extra and marker["pos"] in extra else set()))
        fractions = {major: 1.0 - minor_fraction}
        if minor:
            other = allele_for(marker, set(minor))
            fractions[other] = fractions.get(other, 0.0) + minor_fraction
        evidence[(marker["chrom"], marker["pos"])] = {"fractions": fractions, "depth": depth, "inferred": False}
    return evidence


def classify(evidence, markers=None, families=None):
    default_markers, default_families = barcode()
    return classify_lineage(markers or default_markers, families or default_families, evidence, "pileup", 5)


def test_h37rv_like_isolate_is_lineage_4_9_not_a_fallback():
    result = classify(isolate(["lineage4", "lineage4.9"]))
    assert result["lineage"] == "lineage4.9"
    assert result["lineage_path"] == ["lineage4", "lineage4.9"]
    assert result["confidence"] == "high"
    assert result["matched_snps"] == result["total_snps"] == 10


def test_deepest_consistent_sublineage_is_reported():
    path = ["lineage2", "lineage2.2", "lineage2.2.1", "lineage2.2.1.1"]
    result = classify(isolate(path))
    assert result["lineage"] == "lineage2.2.1.1"
    assert result["lineage_path"] == path
    assert result["family"] == barcode()[1]["lineage2.2.1.1"]
    assert result["confidence"] == "high"
    assert not result["possible_mixed_infection"]


def test_no_coverage_gives_unknown_instead_of_the_first_barcode_lineage():
    result = classify({})
    assert result["lineage"] == "unknown"
    assert result["reason"]


def test_call_does_not_depend_on_barcode_file_order():
    markers, families = barcode()
    shuffled = list(markers)
    random.Random(7).shuffle(shuffled)
    reordered = dict(sorted(families.items(), reverse=True))
    for path in (["lineage4", "lineage4.3", "lineage4.3.4", "lineage4.3.4.2"], ["lineage1", "lineage1.1"]):
        assert classify(isolate(path))["lineage"] == classify(isolate(path), shuffled, reordered)["lineage"] == path[-1]


def test_two_supported_sibling_sublineages_report_their_parent():
    path = ["lineage2", "lineage2.2", "lineage2.2.1", "lineage2.2.1.1", "lineage2.2.1.2"]
    markers, families = barcode()
    for order in (markers, list(reversed(markers))):
        result = classify(isolate(path), order, families)
        assert result["lineage"] == "lineage2.2.1"
        assert result["ambiguous_sublineages"] == ["lineage2.2.1.1", "lineage2.2.1.2"]
        assert result["confidence"] == "low"


def test_a_single_stray_marker_does_not_win():
    markers, _ = barcode()
    stray = next(m["pos"] for m in markers if m["lineage"] == "lineage3.1.1")
    result = classify(isolate(["lineage4", "lineage4.9"], extra={stray}))
    assert result["lineage"] == "lineage4.9"


def test_minor_lineage_is_flagged_as_possible_mixed_infection():
    result = classify(isolate(["lineage2", "lineage2.2"], minor=["lineage4", "lineage4.3"], minor_fraction=0.3))
    assert result["lineage"] == "lineage2.2"
    assert result["possible_mixed_infection"]
    assert {"lineage4", "lineage4.3"} <= {m["lineage"] for m in result["minor_lineages"]}
    assert result["confidence"] != "high"


def test_pileup_reader_keeps_alleles_aligned_with_counts(tmp_path):
    pileup = tmp_path / "s.barcode_pileup.tsv"
    pileup.write_text("NC_000962.3\t1131\tG\tA,<*>\t3,27,0\nNC_000962.3\t4206\tC\t<*>\t40,0\n")
    observed = read_pileup(pileup, {("NC_000962.3", 1131), ("NC_000962.3", 4206)})
    assert observed[("NC_000962.3", 1131)]["fractions"] == pytest.approx({"G": 0.1, "A": 0.9})
    assert observed[("NC_000962.3", 1131)]["depth"] == 30
    assert observed[("NC_000962.3", 4206)]["fractions"] == {"C": 1.0}


def write_vcf(path, records):
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS"]
    for pos, ref, alt in sorted(records):
        lines.append(f"NC_000962.3\t{pos}\t.\t{ref}\t{alt}\t60\tPASS\t.\tGT:AD\t1/1:0,30")
    path.write_text("\n".join(lines) + "\n")


def run_classifier(tmp_path, vcf, pileup=None):
    out = tmp_path / "out.lineage.json"
    args = ["--vcf", vcf, "--barcode", BARCODE, "--reference", REFERENCE, "--sample_id", "S",
            "--output", out, "--min_depth", 5]
    if pileup is not None:
        args += ["--pileup", pileup]
    run_script("lineage_classifier.py", *args)
    return json.loads(out.read_text())


def test_short_vcf_uses_only_the_markers_it_lists(tmp_path):
    markers, _ = barcode()
    ref = reference()["NC_000962.3"]
    records = [(m["pos"], ref[m["pos"] - 1], m["allele"]) for m in markers if m["lineage"] == "lineage2"]
    vcf = tmp_path / "s.vcf"
    write_vcf(vcf, records)
    result = run_classifier(tmp_path, vcf)
    assert result["lineage"] == "lineage2"
    assert result["evidence"] == "vcf"
    assert result["vcf_snps"] == 10
    assert result["confidence"] == "high"
    assert "reference_inferred" not in result


def test_genome_wide_vcf_infers_reference_markers_with_medium_confidence(tmp_path):
    markers, _ = barcode()
    ref = reference()["NC_000962.3"]
    barcode_positions = {m["pos"] for m in markers}
    positions = [p for p in range(100_000, 4_000_000, 50_000) if p not in barcode_positions][:60]
    records = [(p, ref[p - 1], next(b for b in "ACGT" if b != ref[p - 1])) for p in positions]
    vcf = tmp_path / "s.vcf"
    write_vcf(vcf, records)
    result = run_classifier(tmp_path, vcf)
    assert result["lineage"] == "lineage4.9"
    assert result["reference_inferred"] is True
    assert result["confidence"] == "medium"


def test_empty_pileup_is_unknown(tmp_path):
    vcf = tmp_path / "s.vcf"
    write_vcf(vcf, [])
    pileup = tmp_path / "s.barcode_pileup.tsv"
    pileup.write_text("")
    result = run_classifier(tmp_path, vcf, pileup)
    assert result["lineage"] == "unknown"
    assert result["evidence"] == "pileup"
