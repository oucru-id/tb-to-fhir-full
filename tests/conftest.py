"""Shared pytest helpers: run the bin/ scripts on the synthetic fixtures."""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
BIN = REPO / "bin"
FIXTURES = REPO / "tests" / "fixtures"
GOLDEN = REPO / "tests" / "golden"
FIXED_TIMESTAMP = "2026-01-15T08:00:00+00:00"

sys.path.insert(0, str(BIN))
sys.path.insert(0, str(REPO / "tests"))

SAMPLES = sorted(p.name.split(".")[0] for p in (FIXTURES / "vcf").glob("*.annotated_variants.vcf")
                 if not p.name.startswith("LEGACY"))


def script_env():
    env = dict(os.environ)
    env["TBFHIR_FIXED_TIMESTAMP"] = FIXED_TIMESTAMP
    env["PYTHONPATH"] = str(BIN)
    return env


def run_script(name, *args, cwd=None, check=True):
    proc = subprocess.run([sys.executable, str(BIN / name), *map(str, args)],
                          cwd=cwd, env=script_env(), capture_output=True, text=True)
    if check and proc.returncode != 0:
        raise AssertionError(f"{name} exited {proc.returncode}\nstdout:\n{proc.stdout}\n"
                             f"stderr:\n{proc.stderr}")
    return proc


def strip_marker(src, dst):
    lines = Path(src).read_text(encoding="utf-8").splitlines(keepends=True)
    Path(dst).write_text("".join(ln for ln in lines if not ln.startswith("#")), encoding="utf-8")


@pytest.fixture(scope="session")
def staged(tmp_path_factory):
    """All lineage/coverage files staged together, as v1 CREATE_FHIR did."""
    root = tmp_path_factory.mktemp("staged")
    (root / "lineage_data").mkdir()
    (root / "coverage_data").mkdir()
    for f in (FIXTURES / "lineage").glob("*.json"):
        shutil.copy(f, root / "lineage_data" / f.name)
    for f in (FIXTURES / "coverage").glob("*.json"):
        shutil.copy(f, root / "coverage_data" / f.name)
    meta = root / "metadata"
    meta.mkdir()
    for f in (FIXTURES / "metadata").glob("*.csv"):
        strip_marker(f, meta / f.name)
    return root


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))
