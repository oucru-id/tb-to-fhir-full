#!/usr/bin/env python3
"""Regenerate the synthetic test fixtures.

Every value here is invented. Variant rows are copied from the WHO
catalogue annotation table (data/2026-09-07_annotation_table.tsv.gz) so the
annotations look exactly like bcftools output; depths, qualities, names,
NIKs (prefix 99, which is not an Indonesian province code), addresses and
dates are synthetic.

Run from the repository root:  python3 tests/fixtures/make_fixtures.py
"""
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHROM = "NC_000962.3"

# --------------------------------------------------------------------------
# Catalogue rows: (pos, ref, alt, gene, drugs, effect, grades, variant_ids)
# --------------------------------------------------------------------------
RPOB_S450L = (761155, "C", "T", "rpoB", "Rifampicin", "missense_variant",
              "Assoc w R", "rpoB_p.Ser450Leu")
KATG_S315T = (2155168, "C", "G", "katG", "Isoniazid", "missense_variant",
              "Assoc w R", "katG_p.Ser315Thr")
GYRA_D94G = (7582, "A", "G", "gyrA", "Levofloxacin,Moxifloxacin", "missense_variant",
             "Assoc w R,Assoc w R", "gyrA_p.Asp94Gly,gyrA_p.Asp94Gly")
RV0678_LOF = (779014, "C", "T", "Rv0678", "Bedaquiline,Clofazimine",
              "loss_of_function_variant", "Assoc w R,Uncertain significance",
              "Rv0678_LoF,Rv0678_LoF")
RRS_A1401G = (1473246, "A", "G", "rrs", "Amikacin,Capreomycin,Kanamycin,Streptomycin",
              "non_coding_transcript_exon_variant",
              "Assoc w R,Assoc w R,Assoc w R,Uncertain significance",
              "rrs_n.1401A>G,rrs_n.1401A>G,rrs_n.1401A>G,rrs_n.1401A>G")
INHA_C777T = (1673425, "C", "T", "inhA", "Ethionamide,Isoniazid", "upstream_gene_variant",
              "Assoc w R,Assoc w R", "inhA_c.-777C>T,inhA_c.-777C>T")
RPOB_DEL = (761100, "CAATTCATGG", "C", "rpoB", "Rifampicin", "inframe_deletion",
            "Assoc w R", "rpoB_p.Gln432_Asp435delinsHis")
PNCA_LOF = (2288690, "G", "T", "pncA", "Pyrazinamide", "loss_of_function_variant",
            "Assoc w R", "pncA_LoF")
PNCA_INTERIM = (2288681, "T", "C", "pncA", "Pyrazinamide", "stop_lost",
                "Assoc w R - Interim", "pncA_p.Ter187Trpext*?")
EMBB_UNCERTAIN = (4246520, "CA", "AG", "embB", "Ethambutol", "missense_variant",
                  "Uncertain significance", "embB_p.Gln3Arg")
DNAA_NOT_ASSOC = (370, "CC", "TT", "dnaA", "Isoniazid", "missense_variant",
                  "Not assoc w R", "dnaA_p.Pro124Leu")
# The same allele as RV0678_LOF, but as a legacy one-row-per-pair table would
# have left it after bcftools collapsed the grade (two drugs, one grade).
RV0678_COLLAPSED = (779014, "C", "T", "Rv0678", "Bedaquiline,Clofazimine",
                    "loss_of_function_variant", "Assoc w R", "Rv0678_LoF")

UNANNOTATED = (1000000, "G", "A")


def caller_info(depth, style):
    """Caller INFO fields in the three shapes the converter understands."""
    if style == "gatk":      # AF from GATK-style AF
        return f"AC=2;AF=1.00;AN=2;DP={depth};MQ=60.00"
    if style == "dp4":       # bcftools DP4
        ref = depth // 10
        alt = depth - ref
        return f"DP={depth};DP4={ref // 2},{ref - ref // 2},{alt // 2},{alt - alt // 2};MQ=59.12"
    return f"DP={depth}"     # medaka-like: depth only


def record(row, qual=2500.0, depth=80, style="gatk"):
    if len(row) == 3:
        pos, ref, alt = row
        return f"{CHROM}\t{pos}\t.\t{ref}\t{alt}\t{qual:.2f}\tPASS\t{caller_info(depth, style)}"
    pos, ref, alt, gene, drugs, effect, grades, variant_ids = row
    annotation = (f"GENE={gene};DRUG={drugs};EFFECT={effect};WHO_CLASSIFICATION={grades};"
                  f"VARIANT_ID={variant_ids};GENOME_POSITION={pos}")
    return (f"{CHROM}\t{pos}\t.\t{ref}\t{alt}\t{qual:.2f}\tPASS\t"
            f"{caller_info(depth, style)};{annotation}")


VCF_HEADER = """##fileformat=VCFv4.2
##source=TBtoFHIR synthetic test fixture (SYNTHETIC - not patient data)
##contig=<ID=NC_000962.3,length=4411532>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="ref-fwd, ref-rev, alt-fwd, alt-rev">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
##INFO=<ID=DRUG,Number=.,Type=String,Description="Associated drug(s)">
##INFO=<ID=EFFECT,Number=1,Type=String,Description="Predicted effect of mutation">
##INFO=<ID=WHO_CLASSIFICATION,Number=.,Type=String,Description="WHO resistance classification, one per DRUG, same order">
##INFO=<ID=VARIANT_ID,Number=.,Type=String,Description="WHO graded variant identifier, one per DRUG, same order">
##INFO=<ID=GENOME_POSITION,Number=1,Type=Integer,Description="Genome position">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""

# --------------------------------------------------------------------------
# Samples
# --------------------------------------------------------------------------
SAMPLES = {
    # sample_id: (variant records, coverage spec, lineage)
    "RR01": ([record(RPOB_S450L, depth=85), record(UNANNOTATED, qual=812.4, depth=40)],
             "adequate", ("lineage4", "Euro-American")),
    "MDR01": ([record(RPOB_S450L, depth=92), record(KATG_S315T, depth=70, style="dp4")],
              "adequate", ("lineage2", "East-Asian")),
    "PREXDR01": ([record(GYRA_D94G, depth=66), record(RPOB_S450L, depth=81),
                  record(KATG_S315T, depth=77, style="dp4")],
                 "adequate", ("lineage2", "East-Asian")),
    "XDR01": ([record(GYRA_D94G, depth=58), record(RPOB_S450L, depth=63),
               record(RV0678_LOF, qual=901.0, depth=45, style="dp4"),
               record(RRS_A1401G, depth=120), record(KATG_S315T, depth=88)],
              {"inadequate": ["pncA"]}, ("lineage2", "East-Asian")),
    "HR01": ([record(INHA_C777T, depth=74)], "adequate", ("lineage1", "Indo-Oceanic")),
    "INDEL01": ([record(RPOB_DEL, qual=1500.0, depth=55), record(PNCA_LOF, depth=61, style="plain")],
                "adequate", ("lineage4", "Euro-American")),
    "UNC01": ([record(DNAA_NOT_ASSOC, qual=640.0, depth=33),
               record(EMBB_UNCERTAIN, qual=700.0, depth=38)],
              "not_assessed", None),
    "PZA01": ([record(PNCA_INTERIM, depth=52)], {"inadequate": ["rpoB"]},
              ("lineage3", "East-African-Indian")),
    "TB_001": ([record(RPOB_S450L, depth=90)], "adequate", ("lineage4", "Euro-American")),
    "TB_002": ([record(DNAA_NOT_ASSOC, qual=640.0, depth=35)], {"inadequate": ["katG"]},
               ("lineage2", "East-Asian")),
}

# Expected to abort with exit status 2 (collapsed legacy annotation table).
LEGACY = {"LEGACY01": [record(RV0678_COLLAPSED, depth=50)]}

# Gene -> drugs covered by the WHO resistance targets (from make_who_targets.py).
TARGET_GENES = {
    "dnaA": ["Isoniazid"], "gyrB": ["Levofloxacin", "Moxifloxacin"],
    "gyrA": ["Levofloxacin", "Moxifloxacin"], "Rv0010c": ["Isoniazid"],
    "fgd1": ["Clofazimine", "Delamanid"], "mshA": ["Ethionamide", "Isoniazid"],
    "ccsA": ["Amikacin", "Capreomycin", "Kanamycin"], "Rv0565c": ["Ethionamide"],
    "hadA": ["Isoniazid"], "nusG": ["Rifampicin"], "rpoB": ["Rifampicin"],
    "rpoC": ["Rifampicin"], "mmpL5": ["Bedaquiline", "Clofazimine"],
    "mmpS5": ["Bedaquiline", "Clofazimine"], "Rv0678": ["Bedaquiline", "Clofazimine"],
    "rpsL": ["Streptomycin"], "rplC": ["Linezolid"],
    "Rv1129c": ["Isoniazid", "Levofloxacin", "Moxifloxacin", "Rifampicin"],
    "fbiC": ["Clofazimine", "Delamanid"], "sigE": ["Pyrazinamide"],
    "Rv1258c": ["Isoniazid", "Pyrazinamide", "Streptomycin"], "embR": ["Ethambutol"],
    "atpE": ["Bedaquiline"], "rrs": ["Amikacin", "Capreomycin", "Kanamycin", "Streptomycin"],
    "rrl": ["Capreomycin", "Linezolid"], "inhA": ["Ethionamide", "Isoniazid"],
    "rpsA": ["Pyrazinamide"], "tsnR": ["Linezolid"], "tlyA": ["Capreomycin"],
    "bacA": ["Amikacin", "Capreomycin", "Kanamycin", "Streptomycin"],
    "ndh": ["Delamanid", "Ethionamide", "Isoniazid"], "katG": ["Isoniazid"],
    "PPE35": ["Pyrazinamide"], "Rv1979c": ["Bedaquiline", "Clofazimine"],
    "pncA": ["Pyrazinamide"], "eis": ["Amikacin", "Kanamycin"], "ahpC": ["Isoniazid"],
    "Rv2477c": ["Amikacin", "Ethambutol", "Kanamycin", "Levofloxacin", "Moxifloxacin",
                "Rifampicin", "Streptomycin"],
    "pepQ": ["Bedaquiline", "Clofazimine"], "Rv2680": ["Capreomycin"],
    "Rv2681": ["Capreomycin"],
    "Rv2752c": ["Ethambutol", "Isoniazid", "Levofloxacin", "Moxifloxacin", "Rifampicin"],
    "Rv2983": ["Clofazimine", "Delamanid"], "Rv3083": ["Ethionamide"],
    "whiB7": ["Amikacin", "Kanamycin", "Streptomycin"], "Rv3236c": ["Pyrazinamide"],
    "lpqB": ["Bedaquiline", "Rifampicin"], "mtrB": ["Bedaquiline", "Rifampicin"],
    "mtrA": ["Bedaquiline", "Rifampicin"], "fbiA": ["Clofazimine", "Delamanid"],
    "fbiB": ["Clofazimine", "Delamanid"], "rpoA": ["Rifampicin"], "ddn": ["Delamanid"],
    "clpC1": ["Pyrazinamide"], "panD": ["Pyrazinamide"],
    "glpK": ["Ethambutol", "Isoniazid", "Levofloxacin", "Moxifloxacin", "Rifampicin",
             "Streptomycin"],
    "embC": ["Ethambutol"], "embA": ["Ethambutol"], "embB": ["Ethambutol"],
    "aftB": ["Ethambutol"], "ubiA": ["Ethambutol"], "ethA": ["Ethionamide"],
    "ethR": ["Ethionamide"], "whiB6": ["Amikacin", "Capreomycin", "Kanamycin"],
    "gid": ["Streptomycin"],
}


def coverage_summary(sample_id, spec, min_depth=10, min_breadth=0.95):
    """Same shape as scripts/summarize_coverage.py output."""
    if spec == "not_assessed":
        return {
            "sample_id": sample_id, "coverage_assessed": False,
            "reason": "VCF input: no alignment available, so locus coverage could not be assessed",
            "min_depth": min_depth, "min_breadth": min_breadth, "genes": {}, "drugs": {},
        }
    inadequate = set(spec.get("inadequate", [])) if isinstance(spec, dict) else set()
    genes = {}
    for gene in sorted(TARGET_GENES):
        bad = gene in inadequate
        genes[gene] = {
            "bases": 1500,
            "mean_depth": 7.5 if bad else 100.0,
            f"breadth_at_{min_depth}x": 0.42 if bad else 1.0,
            "adequate": not bad,
            "drugs": sorted(TARGET_GENES[gene]),
        }
    drug_to_genes = {}
    for gene, drugs in TARGET_GENES.items():
        for drug in drugs:
            drug_to_genes.setdefault(drug, set()).add(gene)
    drugs = {}
    for drug in sorted(drug_to_genes):
        loci = sorted(drug_to_genes[drug])
        bad = sorted(g for g in loci if g in inadequate)
        drugs[drug] = {"loci": loci, "assessable": not bad, "inadequate_loci": bad,
                       "loci_without_coverage_data": []}
    return {"sample_id": sample_id, "coverage_assessed": True, "min_depth": min_depth,
            "min_breadth": min_breadth, "genes": genes, "drugs": drugs}


def lineage_summary(sample_id, lineage, family):
    return {"lineage": lineage, "confidence": "high", "score": 1.0, "matched_snps": 10,
            "total_snps": 10, "percentage": 100.0, "family": family, "sample_id": sample_id}


# --------------------------------------------------------------------------
# Clinical metadata (synthetic). Rows are present for some samples only, so
# "no clinical row" and "orphan row" cases are covered.
# --------------------------------------------------------------------------
PATIENT_COLUMNS = [
    "sample_id", "family_name", "given_name", "gender", "birth_date", "age", "nik",
    "diagnosis_date", "treatment_history", "contact_number", "address", "city", "state",
    "province_code", "city_code", "district_code", "village_code", "latitude", "longitude",
    "citizenship_status", "specimen_type_code", "specimen_type_display",
    "specimen_collection_method_code", "specimen_collection_method_display",
    "specimen_collection_method_text", "specimen_quantity_value", "specimen_quantity_unit",
    "specimen_collected_date", "specimen_received_date",
]


def patient_row(sample_id, n, **overrides):
    row = {
        "sample_id": sample_id, "family_name": f"Sintetis{n}", "given_name": "Pasien",
        "gender": "laki-laki" if n % 2 else "perempuan", "birth_date": f"19{70 + n}-0{1 + n % 9}-15",
        "age": str(50 - n), "nik": f"99010100000000{n:02d}", "diagnosis_date": "2026-01-05",
        "treatment_history": "new", "contact_number": f"+62-800-0000-{n:04d}",
        "address": f"Jl. Contoh No. {n}, Kec. Uji", "city": "Kota Contoh", "state": "Provinsi Uji",
        "province_code": "99", "city_code": "9901", "district_code": "990101",
        "village_code": "9901011001", "latitude": "-7.0000", "longitude": "110.0000",
        "citizenship_status": "WNI", "specimen_type_code": "119334006",
        "specimen_type_display": "Sputum specimen", "specimen_collection_method_code": "KOFFP",
        "specimen_collection_method_display": "Plate, Cough",
        "specimen_collection_method_text": "Sputum collection", "specimen_quantity_value": "2",
        "specimen_quantity_unit": "mL", "specimen_collected_date": "2026-01-10",
        "specimen_received_date": "2026-01-11",
    }
    row.update(overrides)
    return row


PATIENT_ROWS = [
    patient_row("RR01", 1),
    patient_row("MDR01", 2),
    patient_row("PREXDR01", 3),
    patient_row("XDR01", 4),
    patient_row("HR01", 5),
    # incomplete: no names, no birth date, no specimen dates, no coordinates
    patient_row("INDEL01", 6, family_name="", given_name="", birth_date="",
                specimen_collected_date="", specimen_received_date="",
                latitude="", longitude=""),
    patient_row("TB_001", 7),
    patient_row("TB_002", 8),
    # orphan: no genomic output for this sample
    patient_row("ORPHAN01", 9),
    # UNC01 and PZA01 intentionally have no row
]

ORG_ROW = {
    "org_id": "99000001", "name": "Laboratorium Sintetis", "alias": "LabSin",
    "type_code": "prov", "type_display": "Healthcare Provider", "type_text": "Laboratory",
    "phone": "+62-800-0000-0000", "email": "lab@example.org",
    "address_line": "Jl. Contoh No. 99", "city": "Kota Contoh", "state": "Provinsi Uji",
    "country": "ID", "province_code": "99", "city_code": "9901", "district_code": "990101",
    "latitude": "-7.0001", "longitude": "110.0001",
}

PRACTITIONER_ROW = {
    "practitioner_id": "N99999999", "nik": "9900000000000001", "name": "dr. Praktisi Sintetis",
    "phone": "+62-800-0000-0001", "gender": "female", "birth_date": "1980-01-01",
    "str_kki_number": "SYN-STR-0001", "qualification_period_start": "2020-01-01",
    "role_id": "SYN-ROLE-001", "role_code": "309343006", "role_display": "Physician",
    "role_text": "Clinical microbiologist",
}


# --------------------------------------------------------------------------
# Raw (unannotated) VCFs for the end-to-end `-profile test_vcf` run, as a
# variant caller would write them: GT:DP:GQ per sample, plus lineage barcode
# SNPs so the lineage step has something to call.
# --------------------------------------------------------------------------
REPO = HERE.parents[1]
TEST_DATA = REPO / "tests" / "data"

# Lineages whose barcode alleles differ from H37Rv (the 'lineage4' markers
# equal the reference, so they never appear in a VCF).
RAW_SAMPLES = {
    # sample: (catalogue rows, lineage whose barcode SNPs are added)
    "RR01": ([RPOB_S450L], "lineage1"),
    "MDR01": ([RPOB_S450L, KATG_S315T], "lineage2"),
    "XDR01": ([GYRA_D94G, RPOB_S450L, RV0678_LOF, RRS_A1401G, KATG_S315T], "lineage2"),
    "UNC01": ([DNAA_NOT_ASSOC, EMBB_UNCERTAIN], None),
    "TB_001": ([RPOB_S450L], "lineage1"),
    "TB_002": ([DNAA_NOT_ASSOC], "lineage2"),
}

RAW_HEADER = """##fileformat=VCFv4.2
##source=TBtoFHIR synthetic test input (SYNTHETIC - not patient data)
##contig=<ID=NC_000962.3,length=4411532>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


def read_reference(path):
    seq = []
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)


def lineage_markers(lineage, reference):
    """[(pos, ref, alt)] barcode SNPs of exactly this lineage, from the barcode BED."""
    markers = []
    with open(REPO / "data" / "2026-05-11_lineage.bed") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and fields[3] == lineage:
                pos = int(fields[2])
                assert reference[pos - 1] != fields[4], f"{lineage} marker {pos} equals the reference"
                markers.append((pos, reference[pos - 1], fields[4]))
    return markers


def write_raw_vcfs():
    (TEST_DATA / "vcf").mkdir(parents=True, exist_ok=True)
    reference = read_reference(REPO / "data" / "H37Rv.fasta")
    for sample, (rows, lineage) in RAW_SAMPLES.items():
        records = []
        for pos, ref, alt, *_ann in rows:
            assert reference[pos - 1: pos - 1 + len(ref)] == ref, (sample, pos, ref)
            records.append((pos, ref, alt, 80))
        if lineage:
            records.extend((pos, ref, alt, 60) for pos, ref, alt in lineage_markers(lineage, reference))
        lines = [f"{CHROM}\t{pos}\t.\t{ref}\t{alt}\t1500.00\tPASS\tDP={dp};MQ=60.00\tGT:DP:GQ\t1/1:{dp}:99"
                 for pos, ref, alt, dp in sorted(records)]
        (TEST_DATA / "vcf" / f"{sample}.vcf").write_text(
            RAW_HEADER.format(sample=sample) + "\n".join(lines) + "\n")
    with open(TEST_DATA / "samplesheet_vcf.csv", "w") as handle:
        handle.write("sample,platform,fastq_1,fastq_2,vcf\n")
        for sample in RAW_SAMPLES:
            handle.write(f"{sample},vcf,,,tests/data/vcf/{sample}.vcf\n")


def write_stub_inputs():
    """Placeholder inputs for `-profile test_stub` (stub blocks never read them)."""
    import gzip
    stub = TEST_DATA / "stub"
    stub.mkdir(parents=True, exist_ok=True)
    for name in ("STUBILL_1.fastq.gz", "STUBILL_2.fastq.gz", "STUBONT.fastq.gz"):
        with open(stub / name, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as handle:
            handle.write(b"")
    (stub / "STUBVCF.vcf").write_text(RAW_HEADER.format(sample="STUBVCF"))
    with open(TEST_DATA / "samplesheet_stub.csv", "w") as handle:
        handle.write("sample,platform,fastq_1,fastq_2,vcf\n")
        handle.write("STUBILL,illumina,tests/data/stub/STUBILL_1.fastq.gz,tests/data/stub/STUBILL_2.fastq.gz,\n")
        handle.write("STUBONT,ont,tests/data/stub/STUBONT.fastq.gz,,\n")
        handle.write("STUBVCF,vcf,,,tests/data/stub/STUBVCF.vcf\n")


def write_csv(path, columns, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write("# SYNTHETIC test fixture - invented values, not patient data\n")
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    for sub in ("vcf", "coverage", "lineage", "metadata"):
        (HERE / sub).mkdir(parents=True, exist_ok=True)

    for sample_id, (records, cov_spec, lineage) in {**SAMPLES, **{k: (v, None, None) for k, v in LEGACY.items()}}.items():
        (HERE / "vcf" / f"{sample_id}.annotated_variants.vcf").write_text(
            VCF_HEADER + "\n".join(sorted(records, key=lambda r: int(r.split("\t")[1]))) + "\n")
        if cov_spec is not None:
            (HERE / "coverage" / f"{sample_id}.coverage.json").write_text(
                json.dumps(coverage_summary(sample_id, cov_spec), indent=2) + "\n")
        if lineage:
            (HERE / "lineage" / f"{sample_id}.lineage.json").write_text(
                json.dumps(lineage_summary(sample_id, *lineage), indent=2) + "\n")

    write_csv(HERE / "metadata" / "patient_clinical_metadata.csv", PATIENT_COLUMNS, PATIENT_ROWS)
    write_csv(HERE / "metadata" / "organization_metadata.csv", list(ORG_ROW), [ORG_ROW])
    write_csv(HERE / "metadata" / "practitioner_metadata.csv", list(PRACTITIONER_ROW),
              [PRACTITIONER_ROW])
    write_raw_vcfs()
    write_stub_inputs()
    print(f"fixtures written under {HERE} and {TEST_DATA}")


if __name__ == "__main__":
    main()
