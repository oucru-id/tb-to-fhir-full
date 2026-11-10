#!/usr/bin/env python3
"""Refuse to commit files that look like they carry PHI or credentials.

Checks (each file given on the command line, or every tracked file when
called without arguments):

* CSV files with data rows must start with a ``# SYNTHETIC`` line.
  Header-only templates are fine.
* Spreadsheets (.xlsx/.xls) and SQLite databases are never allowed.
* 16-digit numbers (Indonesian NIK format) are rejected unless they are on
  the synthetic allowlist below or the line says SYNTHETIC.
* Base64-looking blobs of 4 kB or more in source/docs are rejected (an
  embedded report was found in a dev script once).
* JSON with non-empty ``client_secret`` / ``access_token`` / ``refresh_token``
  / ``id_token`` values is rejected.

Exit status 1 lists every offending file and the reason.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

SYNTHETIC_NIKS = {"0000012345678901", "1234567890123456"}
NIK_RE = re.compile(r"(?<![\d.])\d{16}(?![\d.])")
B64_RE = re.compile(r"[A-Za-z0-9+/]{4000,}={0,2}")
SECRET_KEYS = {"client_secret", "access_token", "refresh_token", "id_token"}

TEXT_SUFFIXES = {".py", ".sh", ".nf", ".config", ".md", ".txt", ".json", ".yml",
                 ".yaml", ".csv", ".tsv", ".bed", ".vcf", ".html", ".cff", ""}
B64_SUFFIXES = {".py", ".sh", ".nf", ".config", ".md", ".html"}
FORBIDDEN_SUFFIXES = {".xlsx", ".xls", ".db", ".sqlite", ".sqlite3"}


def tracked_files():
    out = subprocess.run(["git", "ls-files", "-z"], capture_output=True,
                         check=True).stdout.decode()
    return [p for p in out.split("\0") if p]


def check_csv(path, lines):
    content = [ln for ln in lines if ln.strip()]
    if len(content) <= 1:
        return None
    if not content[0].lstrip().startswith("# SYNTHETIC"):
        return "CSV has data rows but no leading '# SYNTHETIC' marker"
    return None


def check_json_secrets(path, text):
    try:
        data = json.loads(text)
    except ValueError:
        return None
    found = []

    def walk(obj):
        if isinstance(obj, dict):
            for key, value in obj.items():
                if key in SECRET_KEYS and isinstance(value, str) and value.strip("x "):
                    found.append(key)
                walk(value)
        elif isinstance(obj, list):
            for item in obj:
                walk(item)

    walk(data)
    return f"JSON carries credential values: {', '.join(sorted(set(found)))}" if found else None


def check_file(name):
    path = Path(name)
    if not path.is_file():
        return []
    suffix = path.suffix.lower()
    if suffix in FORBIDDEN_SUFFIXES:
        return [f"{suffix} files may hold clinical data and are not allowed"]
    if suffix not in TEXT_SUFFIXES:
        return []
    try:
        text = path.read_text(encoding="utf-8-sig")
    except (UnicodeDecodeError, OSError):
        return []
    problems = []
    lines = text.splitlines()
    if suffix == ".csv":
        msg = check_csv(path, lines)
        if msg:
            problems.append(msg)
    if suffix == ".json":
        msg = check_json_secrets(path, text)
        if msg:
            problems.append(msg)
    for lineno, line in enumerate(lines, 1):
        if "SYNTHETIC" in line:
            continue
        for match in NIK_RE.findall(line):
            if match not in SYNTHETIC_NIKS:
                problems.append(f"line {lineno}: 16-digit number (NIK format)")
                break
    if suffix in B64_SUFFIXES and B64_RE.search(text):
        problems.append("contains a base64-like blob of 4 kB or more")
    return problems


def main(argv):
    names = argv[1:] or tracked_files()
    failed = False
    for name in names:
        for problem in check_file(name):
            failed = True
            print(f"{name}: {problem}", file=sys.stderr)
    if failed:
        print("\nRefusing: remove the data, or mark genuinely synthetic fixtures "
              "with a '# SYNTHETIC' first line.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
