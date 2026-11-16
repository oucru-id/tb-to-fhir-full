import json

from conftest import run_script

CHROM = "NC_000962.3"
ANNOTATION = [
    ("761155", "C", "T", "rpoB", "Rifampicin", "Assoc w R", "rpoB_p.Ser450Leu"),
    ("3878506", "G", "A", "rpoA", "Rifampicin", "Uncertain significance", "rpoA_p.Gly31Ser"),
    ("7582", "A", "G", "gyrA", "Levofloxacin,Moxifloxacin", "Assoc w R,Assoc w R - Interim", "gyrA_p.Asp94Gly"),
    ("3086820", "G", "A", "ald", "Cycloserine", "Uncertain significance", "ald_c.-32G>A"),
]


def write_annotation(path):
    lines = []
    for pos, ref, alt, gene, drugs, grades, variant in ANNOTATION:
        n = len(drugs.split(","))
        lines.append("\t".join([CHROM, pos, ref, alt, gene, drugs, "missense_variant", grades,
                                ",".join([variant] * n), pos]))
    path.write_text("\n".join(lines) + "\n")


def targets_by_gene(path):
    rows = {}
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        rows[fields[3]] = fields
    return rows


def test_targets_mark_only_grade_1_2_genes_as_gating(tmp_path):
    annotation = tmp_path / "annotation.tsv"
    write_annotation(annotation)
    bed = tmp_path / "who_targets.bed"
    run_script("make_who_targets.py", "--annotation-table", annotation, "--output", bed)
    rows = targets_by_gene(bed)
    assert bed.read_text().splitlines()[0].split("\t")[-1] == "gating_drugs"
    assert rows["rpoB"][4:] == ["Rifampicin", "Rifampicin"]
    assert rows["rpoA"][4:] == ["Rifampicin", "."]
    assert rows["gyrA"][4:] == ["Levofloxacin;Moxifloxacin", "Levofloxacin;Moxifloxacin"]
    assert rows["ald"][4:] == ["Cycloserine", "."]


def write_coverage(tmp_path, targets, covered):
    regions, thresholds = tmp_path / "s.regions.bed", tmp_path / "s.thresholds.bed"
    region_lines, threshold_lines = [], ["#chrom\tstart\tend\tregion\t10X\t30X"]
    for line in targets.read_text().splitlines():
        if line.startswith("#"):
            continue
        chrom, start, end, gene = line.split("\t")[:4]
        span = int(end) - int(start)
        at_depth = int(span * covered.get(gene, 1.0))
        region_lines.append(f"{chrom}\t{start}\t{end}\t{gene}\t40.0")
        threshold_lines.append(f"{chrom}\t{start}\t{end}\t{gene}\t{at_depth}\t{at_depth}")
    regions.write_text("\n".join(region_lines) + "\n")
    thresholds.write_text("\n".join(threshold_lines) + "\n")
    return regions, thresholds


def summarize(tmp_path, targets, covered):
    regions, thresholds = write_coverage(tmp_path, targets, covered)
    out = tmp_path / "s.coverage.json"
    run_script("summarize_coverage.py", "--sample-id", "S", "--targets", targets, "--regions", regions,
               "--thresholds", thresholds, "--min-depth", 10, "--min-breadth", 0.95, "--output", out)
    return json.loads(out.read_text())


def make_targets(tmp_path):
    annotation = tmp_path / "annotation.tsv"
    write_annotation(annotation)
    bed = tmp_path / "who_targets.bed"
    run_script("make_who_targets.py", "--annotation-table", annotation, "--output", bed)
    return bed


def test_poor_coverage_of_a_non_gating_gene_no_longer_blocks_susceptible(tmp_path):
    summary = summarize(tmp_path, make_targets(tmp_path), {"rpoA": 0.5})
    rifampicin = summary["drugs"]["Rifampicin"]
    assert rifampicin["assessable"] is True
    assert rifampicin["gating_loci"] == ["rpoB"]
    assert rifampicin["inadequate_loci"] == []
    assert rifampicin["inadequate_other_loci"] == ["rpoA"]
    assert summary["gating"] == "loci with WHO grade 1-2 variants"


def test_poor_coverage_of_a_gating_gene_still_gives_indeterminate(tmp_path):
    summary = summarize(tmp_path, make_targets(tmp_path), {"rpoB": 0.9})
    assert summary["drugs"]["Rifampicin"]["assessable"] is False
    assert summary["drugs"]["Rifampicin"]["inadequate_loci"] == ["rpoB"]


def test_drug_without_grade_1_2_variants_is_not_assessable(tmp_path):
    summary = summarize(tmp_path, make_targets(tmp_path), {})
    cycloserine = summary["drugs"]["Cycloserine"]
    assert cycloserine["assessable"] is False
    assert cycloserine["gating_loci"] == []
    assert cycloserine["reason"]


def test_targets_without_gating_column_keep_the_all_loci_rule(tmp_path):
    bed = make_targets(tmp_path)
    old = tmp_path / "old_targets.bed"
    old.write_text("".join("\t".join(line.split("\t")[:5]) + "\n" for line in bed.read_text().splitlines()))
    summary = summarize(tmp_path, old, {"rpoA": 0.5})
    assert summary["drugs"]["Rifampicin"]["assessable"] is False
    assert summary["drugs"]["Rifampicin"]["inadequate_loci"] == ["rpoA"]
    assert summary["gating"] == "all catalogue loci"
