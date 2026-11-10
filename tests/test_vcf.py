import gzip
import shutil

import pytest

from conftest import FIXTURES
from tbfhir_vcf import AnnotationError, genomic_hgvs, read_vcf, resolve_variant_drug_pairs


@pytest.mark.parametrize("pos,ref,alt,expected", [
    (761155, "C", "T", "NC_000962.3:g.761155C>T"),
    (761100, "CAATTCATGG", "C", "NC_000962.3:g.761101_761109del"),
    (100, "AC", "A", "NC_000962.3:g.101del"),
    (761102, "A", "ATTC", "NC_000962.3:g.761102_761103insTTC"),
    (4246520, "CA", "AG", "NC_000962.3:g.4246520_4246521delinsAG"),
    (5, "C", "GA", "NC_000962.3:g.5delinsGA"),
])
def test_genomic_hgvs(pos, ref, alt, expected):
    assert genomic_hgvs("NC_000962.3", pos, ref, alt) == expected


def test_reads_plain_and_gzipped_vcf(tmp_path):
    plain = FIXTURES / "vcf" / "XDR01.annotated_variants.vcf"
    gz = tmp_path / "XDR01.vcf.gz"
    with open(plain, "rb") as src, gzip.open(gz, "wb") as dst:
        shutil.copyfileobj(src, dst)
    a, b = read_vcf(plain), read_vcf(gz)
    assert [(r.POS, r.REF, r.ALT, r.INFO) for r in a] == [(r.POS, r.REF, r.ALT, r.INFO) for r in b]
    assert len(a) == 5


def test_per_drug_grades_stay_aligned():
    record = next(r for r in read_vcf(FIXTURES / "vcf" / "XDR01.annotated_variants.vcf") if r.POS == 779014)
    assert resolve_variant_drug_pairs(record) == [
        ("Bedaquiline", "Assoc w R", "Rv0678_LoF"), ("Clofazimine", "Uncertain significance", "Rv0678_LoF")]


def test_collapsed_legacy_row_is_an_error():
    record = read_vcf(FIXTURES / "vcf" / "LEGACY01.annotated_variants.vcf")[0]
    with pytest.raises(AnnotationError):
        resolve_variant_drug_pairs(record)
    assert resolve_variant_drug_pairs(record, allow_legacy=True) == [
        ("Bedaquiline", "Assoc w R", "Rv0678_LoF"), ("Clofazimine", "Assoc w R", "Rv0678_LoF")]
