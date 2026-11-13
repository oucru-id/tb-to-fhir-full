#!/usr/bin/env python3
"""Generate the v1 equivalence goldens from the baseline commit.

The baseline is the code that was running in production when the 2.0 work
started (origin/main 550a8dc, newer than the v1.6.0 tag). Its FHIR scripts
are extracted with ``git show`` into a temporary directory and run on the
synthetic fixtures exactly the way the v1 pipeline ran them: CREATE_FHIR
received every sample's lineage and coverage file, then MERGE_CLINICAL_FHIR
ran per bundle with all three metadata CSVs.

The baseline scripts need pandas, so pass an interpreter that has it:

    python3 tests/compat/run_v1.py --python ~/medaka-env/bin/python

Outputs (normalised, see tests/normalize.py) go to tests/golden/v1/.
"""
import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests"))
from normalize import dumps, normalize_bundle  # noqa: E402

FIXTURES = REPO / "tests" / "fixtures"
SCRIPTS = ["annotated_to_fhir.py", "merge_clinical_fhir.py", "clinical_metadata_parser.py"]
LEGACY_SAMPLES = {"LEGACY01"}


def baseline_commit():
    return (Path(__file__).with_name("BASELINE").read_text().split()[0]).strip()


def strip_marker(src, dst):
    """v1 reads CSVs with pandas, which cannot skip the '# SYNTHETIC' line."""
    lines = Path(src).read_text(encoding="utf-8").splitlines(keepends=True)
    Path(dst).write_text("".join(ln for ln in lines if not ln.startswith("#")), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--python", default=sys.executable,
                        help="interpreter with pandas for the baseline scripts")
    parser.add_argument("--commit", default=None, help="baseline commit (default: tests/compat/BASELINE)")
    parser.add_argument("--out", default=str(REPO / "tests" / "golden" / "v1"))
    args = parser.parse_args()

    commit = args.commit or baseline_commit()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="tbfhir-v1-") as tmp:
        tmp = Path(tmp)
        scripts = tmp / "scripts"
        scripts.mkdir()
        for name in SCRIPTS:
            blob = subprocess.run(["git", "-C", str(REPO), "show", f"{commit}:scripts/{name}"],
                                  check=True, capture_output=True).stdout
            (scripts / name).write_bytes(blob)

        meta = tmp / "metadata"
        meta.mkdir()
        for name in ("patient_clinical_metadata", "organization_metadata", "practitioner_metadata"):
            strip_marker(FIXTURES / "metadata" / f"{name}.csv", meta / f"{name}.csv")

        # CREATE_FHIR staged every sample's lineage and coverage file together.
        work = tmp / "work"
        (work / "lineage_data").mkdir(parents=True)
        (work / "coverage_data").mkdir()
        for f in (FIXTURES / "lineage").glob("*.lineage.json"):
            shutil.copy(f, work / "lineage_data" / f.name)
        for f in (FIXTURES / "coverage").glob("*.coverage.json"):
            shutil.copy(f, work / "coverage_data" / f.name)

        summary = {}
        for vcf in sorted((FIXTURES / "vcf").glob("*.annotated_variants.vcf")):
            sample = vcf.name.split(".")[0]          # Nextflow simpleName
            staged = work / f"{sample}.vcf"
            shutil.copy(vcf, staged)
            genomic = work / f"{sample}.fhir.json"
            run = subprocess.run(
                [args.python, str(scripts / "annotated_to_fhir.py"),
                 "--input", str(staged), "--output", str(genomic),
                 "--lineage_dir", str(work / "lineage_data"),
                 "--coverage_dir", str(work / "coverage_data"),
                 "--pipeline_version", "1.6.0",
                 "--coverage_min_depth", "10", "--coverage_min_breadth", "0.95",
                 "--organization_metadata", str(meta / "organization_metadata.csv")],
                cwd=work, capture_output=True, text=True)
            summary[sample] = {"genomic_exit": run.returncode}
            if sample in LEGACY_SAMPLES:
                continue
            if run.returncode != 0:
                print(run.stdout, run.stderr, file=sys.stderr)
                sys.exit(f"baseline annotated_to_fhir.py failed for {sample}")

            merged = work / f"{sample}.merged.fhir.json"
            run = subprocess.run(
                [args.python, str(scripts / "merge_clinical_fhir.py"),
                 "--input", str(genomic), "--output", str(merged),
                 "--patient_metadata", str(meta / "patient_clinical_metadata.csv"),
                 "--organization_metadata", str(meta / "organization_metadata.csv"),
                 "--practitioner_metadata", str(meta / "practitioner_metadata.csv"),
                 "--pipeline_version", "1.6.0",
                 "--filter_min_depth", "5", "--filter_min_quality", "20",
                 "--coverage_min_depth", "10", "--coverage_min_breadth", "0.95"],
                cwd=work, capture_output=True, text=True)
            summary[sample]["merge_exit"] = run.returncode
            if run.returncode != 0:
                print(run.stdout, run.stderr, file=sys.stderr)
                sys.exit(f"baseline merge_clinical_fhir.py failed for {sample}")

            for src, suffix in ((genomic, "genomic"), (merged, "merged")):
                bundle = normalize_bundle(json.loads(src.read_text()))
                (out_dir / f"{sample}.{suffix}.json").write_text(dumps(bundle))

        (out_dir / "summary.json").write_text(
            json.dumps({"baseline_commit": commit, "samples": summary}, indent=2, sort_keys=True) + "\n")
    print(f"v1 goldens written to {out_dir} from {commit}")


if __name__ == "__main__":
    main()
