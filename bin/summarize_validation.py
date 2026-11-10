#!/usr/bin/env python3

import argparse
import json
import os
import sys
from collections import defaultdict


def file_of(outcome):
    for ext in outcome.get("extension", []):
        if ext.get("url", "").endswith("operationoutcome-file"):
            return ext.get("valueString")
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--validation", required=True, help="validator -output JSON")
    p.add_argument("--files", nargs="+", required=True, help="the files that were validated")
    p.add_argument("--fail_on_error", action="store_true")
    p.add_argument("--offline", action="store_true",
                   help="validator ran with -tx n/a: errors that only say a code could not be "
                        "checked ('No server available') are counted as warnings")
    p.add_argument("--max_messages", type=int, default=10)
    args = p.parse_args()

    with open(args.validation, encoding="utf-8") as handle:
        data = json.load(handle)
    outcomes = ([e["resource"] for e in data.get("entry", [])]
                if data.get("resourceType") == "Bundle" else [data])

    counts = defaultdict(lambda: {"error": 0, "warning": 0, "information": 0})
    messages = defaultdict(list)
    unchecked = 0
    names = {os.path.basename(f): f for f in args.files}
    for index, outcome in enumerate(outcomes):
        source = file_of(outcome)
        if source is None and len(args.files) == len(outcomes):
            source = args.files[index]
        name = os.path.basename(source or "?")
        for issue in outcome.get("issue", []):
            severity = "error" if issue.get("severity") in ("error", "fatal") else issue.get("severity", "information")
            text = issue.get("details", {}).get("text") or issue.get("diagnostics", "")
            if severity == "error" and args.offline and "No server available" in text:
                severity = "warning"      
                unchecked += 1
            counts[name][severity] = counts[name].get(severity, 0) + 1
            if severity == "error":
                where = (issue.get("expression") or issue.get("location") or ["?"])[0]
                messages[name].append(f"{where}: {text}")

    with open("validation_summary.tsv", "w", encoding="utf-8") as out:
        out.write("file\terrors\twarnings\tinformation\n")
        for name in sorted(names):
            c = counts.get(name, {"error": 0, "warning": 0, "information": 0})
            out.write(f"{name}\t{c['error']}\t{c['warning']}\t{c['information']}\n")
    passed = [names[n] for n in sorted(names) if counts.get(n, {}).get("error", 0) == 0]
    with open("passed.txt", "w", encoding="utf-8") as out:
        out.write("\n".join(os.path.basename(f) for f in passed) + ("\n" if passed else ""))

    failed = sorted(n for n in names if counts.get(n, {}).get("error", 0))
    print(f"validated {len(names)} file(s): {len(passed)} without errors, {len(failed)} with errors")
    if unchecked:
        print(f"note: {unchecked} code check(s) need a terminology server and were not performed "
              f"(offline validation, -tx n/a)")
    if failed:
        for name in failed:
            print(f"  {name}: {counts[name]['error']} error(s)", file=sys.stderr)
            for msg in messages[name][: args.max_messages]:
                print(f"    - {msg}", file=sys.stderr)
        if args.fail_on_error:
            sys.exit("ERROR: FHIR validation failed (set --fail_on_validation_error false to continue)")


if __name__ == "__main__":
    main()
