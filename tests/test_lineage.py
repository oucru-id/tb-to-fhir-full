"""Lineage calling from barcode SNPs (classifier unchanged from v1.6)."""
from conftest import REPO
from lineage_classifier import classify_lineage, load_barcode_bed

BARCODE = REPO / "data" / "2026-05-11_lineage.bed"


def sample_with(lineage_snps, lineage):
    return {snp["snp_key"]: {} for snp in lineage_snps[lineage]}


def test_lineage_with_full_barcode_support():
    lineage_snps, info = load_barcode_bed(str(BARCODE))
    result = classify_lineage(sample_with(lineage_snps, "lineage2"), lineage_snps, info)
    assert result["lineage"] == "lineage2"
    assert result["confidence"] == "high"
    assert result["matched_snps"] == result["total_snps"] == 10


def test_no_barcode_snp_keeps_v1_behaviour():
    """Documented limitation (README): with no matched barcode SNP, v1.6 -- and 2.0 --
    report the first lineage of the barcode file with score 0 and confidence 'low'."""
    lineage_snps, info = load_barcode_bed(str(BARCODE))
    result = classify_lineage({}, lineage_snps, info)
    assert result["lineage"] == next(iter(lineage_snps))
    assert result["matched_snps"] == 0
    assert result["confidence"] == "low"
