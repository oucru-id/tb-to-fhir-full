#!/usr/bin/env python3

import argparse
import copy
import glob
import json
import os
import sys

from tbfhir_common import bundle_id, full_url, now_iso, today_compact, walk_references
from tbfhir_profiles import get_profile

EXTERNAL_TYPES = {"Location"}
ORDER = ["Device", "Organization", "Practitioner", "PractitionerRole", "Patient", "Specimen",
         "Encounter", "ServiceRequest", "DiagnosticReport", "Observation", "Provenance"]


class LinkError(Exception):
    pass


def load_bundle(path):
    with open(path, encoding="utf-8") as handle:
        bundle = json.load(handle)
    if bundle.get("resourceType") != "Bundle":
        raise LinkError(f"{path}: not a FHIR Bundle")
    return [e["resource"] for e in bundle.get("entry", [])]


def by_type(resources, rtype):
    return [r for r in resources if r["resourceType"] == rtype]


def ref(resource):
    return f"{resource['resourceType']}/{resource['id']}"


def link_sample(genomic, clinical, shared, external, profile, request_method, timestamp):
    genomic = copy.deepcopy(genomic)
    clinical = copy.deepcopy(clinical)
    shared = copy.deepcopy(shared)

    reports = by_type(genomic, "DiagnosticReport")
    if len(reports) != 1:
        raise LinkError(f"genomic bundle has {len(reports)} DiagnosticReports, expected 1")
    report = reports[0]
    fid = report["id"][: -len("-genomic-report")]
    org = next(iter(by_type(shared, "Organization")), None)
    practitioner = next(iter(by_type(shared, "Practitioner")), None)
    service_request = next(iter(by_type(clinical, "ServiceRequest")), None)
    encounter = next(iter(by_type(clinical, "Encounter")), None)
    sample_external = external.get(fid, {})

    if service_request:
        report["basedOn"] = [{"reference": ref(service_request)}]
    if encounter or "Encounter" in sample_external:
        report["encounter"] = {"reference": f"Encounter/{fid}-encounter",
                               "display": "TB Testing Encounter"}
    if org:
        org_ref = {"reference": ref(org)}
        if org.get("name"):
            org_ref["display"] = org["name"]
        performers = report.setdefault("performer", [])
        existing = next((p for p in performers if p.get("reference") == ref(org)), None)
        if existing is None:
            performers.insert(0, org_ref)
        elif "display" not in existing and "display" in org_ref:
            existing["display"] = org_ref["display"]
        if "identifier" not in report:
            identifier = profile.identifier("diagnostic-report", f"TB-GEN-{fid}-{today_compact()}",
                                            org["id"])
            if identifier:
                report["identifier"] = [identifier]
        for obs in by_type(genomic, "Observation"):
            obs.setdefault("performer", [{"reference": ref(org)}])
        for prov in by_type(genomic, "Provenance"):
            agents = prov.setdefault("agent", [])
            agent = next((a for a in agents if a.get("who", {}).get("reference") == ref(org)), None)
            if agent is None:
                agents.append({"type": {"coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/provenance-participant-type",
                    "code": "performer", "display": "Performer"}]}, "who": dict(org_ref)})
            elif "display" not in agent["who"] and "display" in org_ref:
                agent["who"]["display"] = org_ref["display"]
    if practitioner:
        p_ref = {"reference": ref(practitioner)}
        name = (practitioner.get("name") or [{}])[0].get("text")
        if name:
            p_ref["display"] = name
        if not any(p.get("reference") == ref(practitioner) for p in report.get("performer", [])):
            report.setdefault("performer", []).append(p_ref)
    if "performer" in report and not report["performer"]:
        del report["performer"]

    resources = shared + clinical + genomic
    resources = [r for r in resources
                 if not (r["resourceType"] in sample_external and r["id"].startswith(fid))]
    resources.sort(key=lambda r: ORDER.index(r["resourceType"]) if r["resourceType"] in ORDER else len(ORDER))

    index = {ref(r): full_url(r["resourceType"], r["id"]) for r in resources}
    full_urls = set(index.values())
    unresolved = []

    def resolve(value):
        if value in index:
            return index[value]
        if value in full_urls:          
            return value
        rtype, _, rid = value.partition("/")
        if rtype in sample_external and rid.startswith(fid):
            return f"{rtype}/{sample_external[rtype]}"
        if rtype in EXTERNAL_TYPES or value.startswith(("http://", "https://", "#")):
            return value
        unresolved.append(value)
        return value

    for resource in resources:
        walk_references(resource, resolve)
    if unresolved:
        raise LinkError(f"sample {fid}: unresolved reference(s): {', '.join(sorted(set(unresolved)))}")

    entries = []
    for resource in resources:
        entry = {"fullUrl": index[ref(resource)], "resource": resource}
        if request_method == "POST":
            request = {"method": "POST", "url": resource["resourceType"]}
            nik = next((i["value"] for i in resource.get("identifier", [])
                        if i.get("system", "").endswith("/id/nik")), None)
            if nik and resource["resourceType"] in ("Patient", "Practitioner"):
                request["ifNoneExist"] = f"identifier=https://fhir.kemkes.go.id/id/nik|{nik}"
            resource.pop("id", None)
        else:
            request = {"method": "PUT", "url": ref(resource)}
        entry["request"] = request
        entries.append(entry)

    return fid, {
        "resourceType": "Bundle",
        "id": bundle_id("transaction", fid),
        "meta": profile.bundle_meta(timestamp),
        "type": "transaction",
        "timestamp": timestamp,
        "entry": entries,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--genomic', nargs='+', required=True, help='genomic bundle JSON file(s)')
    p.add_argument('--clinical_dir', required=True, help='output directory of clinical_to_fhir.py')
    p.add_argument('--fhir_profile', default='generic', choices=['generic', 'satusehat'])
    p.add_argument('--identifier_system_base', default=None)
    p.add_argument('--external_refs', help='JSON {fhir_id: {"Patient": "<server id>", ...}}')
    p.add_argument('--bundle_request_method', default='PUT', choices=['PUT', 'POST'])
    p.add_argument('--require_clinical', action='store_true',
                   help='fail when a sample has no clinical row')
    p.add_argument('--outdir', default='.')
    args = p.parse_args(argv)

    profile = get_profile(args.fhir_profile, args.identifier_system_base)
    shared_path = os.path.join(args.clinical_dir, "shared.clinical.json")
    shared = load_bundle(shared_path) if os.path.exists(shared_path) else []
    external = {}
    for path in (os.path.join(args.clinical_dir, "external_refs.json"), args.external_refs):
        if path and os.path.exists(path):
            with open(path, encoding="utf-8") as handle:
                for fid, refs in json.load(handle).items():
                    external.setdefault(fid, {}).update(refs)

    clinical = {}
    for path in glob.glob(os.path.join(args.clinical_dir, "*.clinical.json")):
        if os.path.basename(path) == "shared.clinical.json":
            continue
        resources = load_bundle(path)
        patient = next(iter(by_type(resources, "Patient")), None)
        if patient:
            clinical[patient["id"][: -len("-patient")]] = resources

    os.makedirs(args.outdir, exist_ok=True)
    timestamp = now_iso()
    summary = []
    errors = []
    for path in sorted(args.genomic):
        genomic = load_bundle(path)
        report = next(iter(by_type(genomic, "DiagnosticReport")), None)
        if report is None:
            errors.append(f"{path}: no DiagnosticReport")
            continue
        fid = report["id"][: -len("-genomic-report")]
        sample_clinical = clinical.get(fid)
        if sample_clinical is None and "Patient" not in external.get(fid, {}):
            summary.append((fid, "missing_clinical", ""))
            continue
        try:
            fid, bundle = link_sample(genomic, sample_clinical or [], shared, external, profile,
                                      args.bundle_request_method, timestamp)
        except LinkError as exc:
            errors.append(str(exc))
            summary.append((fid, "error", ""))
            continue
        out = os.path.join(args.outdir, f"{fid}.bundle.fhir.json")
        with open(out, "w", encoding="utf-8") as handle:
            json.dump(bundle, handle, indent=2, ensure_ascii=False)
        summary.append((fid, "linked", os.path.basename(out)))

    with open(os.path.join(args.outdir, "link_summary.tsv"), "w", encoding="utf-8") as handle:
        handle.write("fhir_id\tstatus\tbundle\n")
        for row in summary:
            handle.write("\t".join(row) + "\n")

    missing = [r[0] for r in summary if r[1] == "missing_clinical"]
    linked = sum(1 for r in summary if r[1] == "linked")
    print(f"linked {linked} sample(s); {len(missing)} without clinical data")
    if missing:
        print(f"WARNING: no clinical data, genomic bundle only: {', '.join(missing)}")
    if errors:
        sys.exit("ERROR:\n  " + "\n  ".join(errors))
    if missing and args.require_clinical:
        sys.exit(f"ERROR: --require_clinical: no clinical data for {', '.join(missing)}")


if __name__ == "__main__":
    main()
