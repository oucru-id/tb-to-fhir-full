#!/usr/bin/env python3

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


def http_post(url, body, headers, timeout):
    request = urllib.request.Request(url, data=body, headers=headers, method="POST")
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return response.status, response.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read().decode("utf-8", "replace")
    except urllib.error.URLError as exc:
        sys.exit(f"ERROR: cannot reach {url}: {exc.reason}")


def fetch_access_token(auth_base_url, client_id, client_secret):
    token_url = f"{auth_base_url.rstrip('/')}/protocol/openid-connect/token"
    form = urllib.parse.urlencode({"grant_type": "client_credentials", "client_id": client_id,
                                   "client_secret": client_secret, "scope": "openid"}).encode()
    status, body = http_post(token_url, form, {"Content-Type": "application/x-www-form-urlencoded"}, 30)
    if status != 200:
        sys.exit(f"ERROR: token request failed (HTTP {status})")
    try:
        token = json.loads(body).get("access_token")
    except ValueError:
        token = None
    if not token:
        sys.exit("ERROR: token response has no access_token")
    return token


def token_from_file(path):
    path = Path(path)
    if not path.exists():
        sys.exit(f"ERROR: token file not found: {path}")
    data = json.loads(path.read_text())
    token = data.get("access_token")
    if not token:
        sys.exit(f"ERROR: no access_token in {path}")
    expires_at = data.get("expires_at")
    if expires_at:
        expiry = datetime.fromisoformat(expires_at)
        if datetime.now(expiry.tzinfo) >= expiry:
            sys.exit(f"ERROR: token in {path} expired at {expires_at}; run get_access_token.py again")
    return token


def summarise_response(body):
    try:
        data = json.loads(body)
    except (ValueError, TypeError):
        return {"unparsed": True}
    if data.get("resourceType") == "Bundle":
        return {"entries": [{"status": e.get("response", {}).get("status"),
                             "location": e.get("response", {}).get("location")}
                            for e in data.get("entry", [])]}
    if data.get("resourceType") == "OperationOutcome":
        return {"issues": [{"severity": i.get("severity"), "code": i.get("code")}
                           for i in data.get("issue", [])]}
    return {"resourceType": data.get("resourceType")}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fhir_file", required=True, help="FHIR Bundle JSON to upload")
    parser.add_argument("--fhir_server_url", required=True, help="FHIR server base URL (POST target)")
    parser.add_argument("--auth_base_url", default="", help="OAuth/OIDC realm URL for client credentials")
    parser.add_argument("--token_file", default="", help="access_token.json from get_access_token.py")
    parser.add_argument("--max_retries", type=int, default=5)
    args = parser.parse_args()

    if not args.fhir_server_url:
        sys.exit("ERROR: --fhir_server_url is empty")

    client_id = os.environ.get("TBFHIR_CLIENT_ID", "")
    client_secret = os.environ.get("TBFHIR_CLIENT_SECRET", "")
    static_token = os.environ.get("TBFHIR_STATIC_TOKEN", "")
    api_key = os.environ.get("TBFHIR_API_KEY", "")

    if args.token_file:
        token = token_from_file(args.token_file)
    elif static_token:
        token = static_token
    elif args.auth_base_url and client_id and client_secret:
        token = fetch_access_token(args.auth_base_url, client_id, client_secret)
    else:
        sys.exit("ERROR: no credentials: set TBFHIR_CLIENT_ID and TBFHIR_CLIENT_SECRET (with "
                 "--auth_base_url), or TBFHIR_STATIC_TOKEN, or pass --token_file")

    body = Path(args.fhir_file).read_bytes()
    json.loads(body)
    headers = {"Content-Type": "application/fhir+json", "Accept": "application/fhir+json",
               "Authorization": f"Bearer {token}"}
    if api_key:
        headers["X-Api-Key"] = api_key

    print(f"Uploading {os.path.basename(args.fhir_file)} to {args.fhir_server_url}")
    for attempt in range(1, args.max_retries + 1):
        status, response_body = http_post(args.fhir_server_url, body, headers, 120)
        if status != 429:
            break
        wait = 2 ** attempt
        print(f"  HTTP 429 (attempt {attempt}/{args.max_retries}); retrying in {wait}s")
        time.sleep(wait)

    success = 200 <= status < 300
    result = {
        "status": "success" if success else "failed",
        "http_status": status,
        "file": os.path.basename(args.fhir_file),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "server": args.fhir_server_url,
        "response": summarise_response(response_body),
    }
    output = Path(args.fhir_file).name.replace(".fhir.json", "").replace(".json", "") + ".upload.json"
    with open(output, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
    print(f"Upload finished with HTTP {status}")
    if not success:
        sys.exit(f"ERROR: upload failed (HTTP {status}); see {output}")


if __name__ == "__main__":
    main()
