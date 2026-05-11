#!/usr/bin/env python3
"""
  Single NIK lookup:
  python scripts/get_patient_by_nik.py --nik 1234

  Multiple NIKs:
  python scripts/get_patient_by_nik.py --nik 1234 5678

  Output results to a JSON file:
  python scripts/get_patient_by_nik.py --nik 1234 --out results/patient_lookup.json
"""
import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

BASEDIR = Path(__file__).parent.parent
DEFAULT_TOKEN_FILE = str(BASEDIR / "data" / "access_token.json")
DEFAULT_FHIR_BASE_URL = ""
NIK_SYSTEM = "https://fhir.kemkes.go.id/id/nik"


def load_access_token(token_file: str) -> str:
    token_path = Path(token_file)
    if not token_path.exists():
        print(f"ERROR: Token file not found: {token_path}", file=sys.stderr)
        sys.exit(1)

    with open(token_path, "r") as fh:
        token_data = json.load(fh)

    access_token = token_data.get("access_token")
    if not access_token:
        print(f"ERROR: No 'access_token' field found in {token_path}", file=sys.stderr)
        sys.exit(1)

    expires_at = token_data.get("expires_at")
    if expires_at:
        try:
            expiry = datetime.fromisoformat(expires_at)
            if datetime.now(expiry.tzinfo) >= expiry:
                print(
                    f"WARNING: Token in {token_path} expired at {expires_at}. "
                    "Re-run get_access_token.py to refresh it.",
                    file=sys.stderr,
                )
        except ValueError:
            pass

    return access_token


def search_patient_by_nik(fhir_base_url: str, nik: str, access_token: str, timeout: int) -> dict:

    params = urlencode({"identifier": f"{NIK_SYSTEM}|{nik}"})
    url = f"{fhir_base_url.rstrip('/')}/Patient?{params}"

    request = Request(
        url,
        headers={
            "Accept": "application/fhir+json",
            "Authorization": f"Bearer {access_token}",
        },
        method="GET",
    )

    try:
        with urlopen(request, timeout=timeout) as response:
            body = response.read().decode("utf-8", errors="replace")
            status_code = response.getcode()
    except HTTPError as exc:
        error_body = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"FHIR request failed ({exc.code}): {error_body}") from exc
    except URLError as exc:
        raise RuntimeError(f"FHIR request failed: {exc.reason}") from exc

    if status_code < 200 or status_code >= 300:
        raise RuntimeError(f"FHIR request failed ({status_code}): {body}")

    try:
        return json.loads(body)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"FHIR response is not valid JSON: {body}") from exc


def extract_patients(bundle: dict) -> list[dict]:
    
    if bundle.get("resourceType") != "Bundle":
        raise RuntimeError(
            f"Unexpected resourceType '{bundle.get('resourceType')}'; expected Bundle."
        )

    entries = bundle.get("entry", [])
    patients = []
    for entry in entries:
        resource = entry.get("resource", {})
        if resource.get("resourceType") != "Patient":
            continue

        uuid = resource.get("id", "")

        nik_value = ""
        for identifier in resource.get("identifier", []):
            if identifier.get("system") == NIK_SYSTEM:
                nik_value = identifier.get("value", "")
                break

        patients.append({"uuid": uuid, "nik": nik_value})

    return patients


def main():
    parser = argparse.ArgumentParser(
        description="Query a FHIR server to find Patient UUID(s) by NIK identifier"
    )
    parser.add_argument(
        "--nik",
        nargs="+",
        required=True,
        help="One or more NIK values to look up",
    )
    parser.add_argument(
        "--fhir_base_url",
        default=DEFAULT_FHIR_BASE_URL,
        help=f"FHIR server base URL (default: {DEFAULT_FHIR_BASE_URL})",
    )
    parser.add_argument(
        "--token_file",
        default=DEFAULT_TOKEN_FILE,
        help=f"Path to access_token.json (default: {DEFAULT_TOKEN_FILE})",
    )
    parser.add_argument(
        "--static_token",
        default="",
        help="Pre-existing Bearer token string (skips reading token_file)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=30,
        help="Request timeout in seconds (default: 30)",
    )
    parser.add_argument(
        "--out",
        default="",
        help="Optional path to write results as JSON",
    )
    args = parser.parse_args()

    if args.static_token:
        access_token = args.static_token
    else:
        access_token = load_access_token(args.token_file)

    all_results = []

    for nik in args.nik:
        print(f"Searching for NIK: {nik} ...")
        try:
            bundle = search_patient_by_nik(args.fhir_base_url, nik, access_token, args.timeout)
        except RuntimeError as exc:
            print(f"  ERROR: {exc}", file=sys.stderr)
            all_results.append({"nik": nik, "error": str(exc)})
            continue

        try:
            patients = extract_patients(bundle)
        except RuntimeError as exc:
            print(f"  ERROR: {exc}", file=sys.stderr)
            all_results.append({"nik": nik, "error": str(exc)})
            continue

        total = bundle.get("total", len(patients))

        if not patients:
            print(f"  No patient found for NIK {nik} (server total={total})")
            all_results.append({"nik": nik, "total": total, "patients": []})
        else:
            for p in patients:
                print(f"  UUID  : {p['uuid']}")
                print(f"  NIK   : {p['nik']}")
            if len(patients) > 1:
                print(f"  WARNING: {len(patients)} patients matched NIK {nik}. Expected 1.")
            all_results.append({"nik": nik, "total": total, "patients": patients})

    if args.out:
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as fh:
            json.dump(all_results, fh, indent=2)
        print(f"\nResults written to {out_path}")

    errors = [r for r in all_results if "error" in r or not r.get("patients")]
    if errors:
        sys.exit(1)


if __name__ == "__main__":
    main()
