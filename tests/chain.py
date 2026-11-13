"""Run the 2.0 chain (genomic -> clinical -> link) on the synthetic fixtures."""
import json
import shutil
from pathlib import Path

from conftest import FIXTURES, run_script

ORG_ID = "99000001"
DEVICE_ARGS = ["--pipeline_version", "1.6.0", "--filter_min_depth", 5, "--filter_min_quality", 20,
               "--coverage_min_depth", 10, "--coverage_min_breadth", 0.95]


def run_genomic(sample, outdir, profile="generic", performer=True, extra=()):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    vcf = outdir / f"{sample}.annotated_variants.vcf"
    shutil.copy(FIXTURES / "vcf" / vcf.name, vcf)
    args = ["--input", vcf, "--sample_id", sample, "--output", outdir / f"{sample}.genomic.fhir.json",
            "--fhir_profile", profile, *DEVICE_ARGS, *extra]
    lineage = FIXTURES / "lineage" / f"{sample}.lineage.json"
    coverage = FIXTURES / "coverage" / f"{sample}.coverage.json"
    if lineage.exists():
        args += ["--lineage_json", lineage]
    if coverage.exists():
        args += ["--coverage_json", coverage]
    if performer:
        args += ["--performer_org_id", ORG_ID]
    run_script("tbfhir_genomic.py", *args, cwd=outdir)
    return outdir / f"{sample}.genomic.fhir.json"


def run_clinical(samples, outdir, profile="generic", organization=True, practitioner=True, check=True):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ids = outdir / "sample_ids.txt"
    ids.write_text("\n".join(samples) + "\n")
    args = ["--patient_metadata", FIXTURES / "metadata" / "patient_clinical_metadata.csv",
            "--sample_ids", ids, "--fhir_profile", profile, "--outdir", outdir]
    if organization:
        args += ["--organization_metadata", FIXTURES / "metadata" / "organization_metadata.csv"]
    if practitioner:
        args += ["--practitioner_metadata", FIXTURES / "metadata" / "practitioner_metadata.csv"]
    return run_script("clinical_to_fhir.py", *args, cwd=outdir, check=check)


def run_link(genomic_files, clinical_dir, outdir, profile="generic", extra=(), check=True):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    return run_script("link_fhir_bundle.py", "--genomic", *genomic_files, "--clinical_dir", clinical_dir,
                      "--fhir_profile", profile, "--outdir", outdir, *extra, cwd=outdir, check=check)


def run_all(samples, root, profile="generic"):
    root = Path(root)
    genomic = [run_genomic(s, root / "genomic", profile) for s in samples]
    run_clinical(samples, root / "clinical", profile)
    run_link(genomic, root / "clinical", root / "bundle", profile)
    return root


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))
