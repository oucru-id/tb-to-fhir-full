import gzip
import json
import shutil

from conftest import FIXTURES, run_script

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


POSITION_ANNOTATION = [
    ("761155", "C", "T", "rpoB", "Rifampicin", "missense_variant", "Assoc w R", "rpoB_p.Ser450Leu"),
    ("1473246", "A", "G", "rrs", "Amikacin,Kanamycin", "non_coding_transcript_exon_variant",
     "Assoc w R,Assoc w R", "rrs_n.1401A>G"),
    ("1473329", "G", "T", "rrs", "Amikacin,Kanamycin", "non_coding_transcript_exon_variant",
     "Assoc w R,Assoc w R", "rrs_n.1484G>T"),
    ("2155168", "C", "G", "katG", "Isoniazid", "missense_variant", "Assoc w R", "katG_p.Ser315Thr"),
    ("2155500", "GC", "G", "katG", "Isoniazid", "frameshift_variant", "Assoc w R", "katG_p.Ala205fs"),
    ("2153900", "T", "C", "katG", "Isoniazid", "missense_variant", "Uncertain significance", "katG_p.Lys731Glu"),
]


def position_targets(tmp_path):
    annotation = tmp_path / "annotation.tsv"
    annotation.write_text("".join("\t".join([CHROM, pos, ref, alt, gene, drugs, effect, grades,
                                              ",".join([variant] * len(drugs.split(","))), pos]) + "\n"
                                  for pos, ref, alt, gene, drugs, effect, grades, variant in POSITION_ANNOTATION))
    targets, positions = tmp_path / "who_targets.bed", tmp_path / "who_positions.bed"
    run_script("make_who_targets.py", "--annotation-table", annotation, "--output", targets,
               "--positions-output", positions)
    return targets, positions


def position_rows(path):
    return [line.split("\t") for line in path.read_text().splitlines() if not line.startswith("#")]


def summarize_with_positions(tmp_path, targets, positions, gene_breadth, uncovered=()):
    regions, thresholds = write_coverage(tmp_path, targets, gene_breadth)
    position_thresholds = tmp_path / "s.positions.thresholds.bed.gz"
    lines = ["#chrom\tstart\tend\tregion\t10X"]
    for chrom, start, end, gene, _drugs in position_rows(positions):
        span = int(end) - int(start)
        lines.append(f"{chrom}\t{start}\t{end}\t{gene}\t{0 if int(end) in uncovered else span}")
    with gzip.open(position_thresholds, "wt") as handle:
        handle.write("\n".join(lines) + "\n")
    out = tmp_path / "s.coverage.json"
    run_script("summarize_coverage.py", "--sample-id", "S", "--targets", targets, "--regions", regions,
               "--thresholds", thresholds, "--positions", positions, "--position-thresholds",
               position_thresholds, "--min-depth", 10, "--min-breadth", 0.95, "--output", out)
    return json.loads(out.read_text())


def test_positions_bed_lists_graded_positions_of_loci_without_loss_of_function_variants(tmp_path):
    _, positions = position_targets(tmp_path)
    rows = position_rows(positions)
    assert [r[1:] for r in rows] == [["761154", "761155", "rpoB", "Rifampicin"],
                                     ["1473245", "1473246", "rrs", "Amikacin;Kanamycin"],
                                     ["1473328", "1473329", "rrs", "Amikacin;Kanamycin"]]


def test_covered_graded_positions_make_a_partly_covered_locus_assessable(tmp_path):
    targets, positions = position_targets(tmp_path)
    summary = summarize_with_positions(tmp_path, targets, positions, {"rrs": 0.8})
    amikacin = summary["drugs"]["Amikacin"]
    assert amikacin["assessable"] is True
    assert amikacin["position_gated_loci"] == ["rrs"]
    assert summary["genes"]["rrs"]["adequate"] is False
    assert summary["genes"]["rrs"]["graded_positions_at_min_depth"] == 2


def test_uncovered_graded_position_gives_indeterminate_and_is_named(tmp_path):
    targets, positions = position_targets(tmp_path)
    summary = summarize_with_positions(tmp_path, targets, positions, {}, uncovered={1473246})
    for drug in ("Amikacin", "Kanamycin"):
        assert summary["drugs"][drug]["assessable"] is False
        assert summary["drugs"][drug]["inadequate_loci"] == ["rrs"]
        assert summary["drugs"][drug]["uncovered_positions"] == {"rrs": ["1473246"]}


def test_loss_of_function_locus_still_needs_the_whole_gene(tmp_path):
    targets, positions = position_targets(tmp_path)
    summary = summarize_with_positions(tmp_path, targets, positions, {"katG": 0.8})
    isoniazid = summary["drugs"]["Isoniazid"]
    assert isoniazid["assessable"] is False
    assert isoniazid["inadequate_loci"] == ["katG"]
    assert isoniazid["position_gated_loci"] == []


def test_indeterminate_text_names_the_uncovered_graded_position(tmp_path):
    coverage = json.loads((FIXTURES / "coverage" / "TB_002.coverage.json").read_text())
    coverage["drugs"]["Amikacin"].update(assessable=False, inadequate_loci=["rrs"],
                                         uncovered_positions={"rrs": ["1473246"]})
    coverage_path = tmp_path / "TB_002.coverage.json"
    coverage_path.write_text(json.dumps(coverage))
    vcf = tmp_path / "TB_002.annotated_variants.vcf"
    shutil.copy(FIXTURES / "vcf" / vcf.name, vcf)
    out = tmp_path / "TB_002.genomic.fhir.json"
    run_script("tbfhir_genomic.py", "--input", vcf, "--sample_id", "TB_002", "--output", out,
               "--coverage_json", coverage_path, cwd=tmp_path)
    bundle = json.loads(out.read_text())
    panel = next(e["resource"] for e in bundle["entry"]
                 if any(c.get("code") == "89486-5" for c in e["resource"].get("code", {}).get("coding", [])))
    amikacin = next(c for c in panel["component"] if c["code"]["coding"][0]["code"] == "89484-0")
    assert "insufficient coverage at: rrs (graded position 1473246)" in amikacin["valueCodeableConcept"]["text"]


def test_targets_without_gating_column_keep_the_all_loci_rule(tmp_path):
    bed = make_targets(tmp_path)
    old = tmp_path / "old_targets.bed"
    old.write_text("".join("\t".join(line.split("\t")[:5]) + "\n" for line in bed.read_text().splitlines()))
    summary = summarize(tmp_path, old, {"rpoA": 0.5})
    assert summary["drugs"]["Rifampicin"]["assessable"] is False
    assert summary["drugs"]["Rifampicin"]["inadequate_loci"] == ["rpoA"]
    assert summary["gating"] == "all catalogue loci"
