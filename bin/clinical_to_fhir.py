#!/usr/bin/env python3

import argparse
import json
import os
import re
import sys

from clinical_metadata_parser import (MetadataError, get_clinical_value, load_clinical_metadata,
                                      load_organization_metadata, load_practitioner_metadata)
from tbfhir_common import bundle_id, fhir_id, full_url, lookup_sample, now_iso
from tbfhir_profiles import KEMKES_NIK, KEMKES_SD, KEMKES_SYS_IDS, get_profile

FHIR_DATE = re.compile(r"^\d{4}(-(0[1-9]|1[0-2])(-(0[1-9]|[12]\d|3[01]))?)?$")
FHIR_DATETIME = re.compile(
    r"^\d{4}-(0[1-9]|1[0-2])-(0[1-9]|[12]\d|3[01])T([01]\d|2[0-3]):[0-5]\d:([0-5]\d|60)"
    r"(\.\d+)?(Z|[+-]((0\d|1[0-3]):[0-5]\d|14:00))$")
NIK = re.compile(r"^\d{16}$")


def fhir_date(value):
    if not value:
        return None
    value = value.strip()
    if " " in value and FHIR_DATE.match(value.split(" ", 1)[0]):
        value = value.split(" ", 1)[0]
    return value if FHIR_DATE.match(value) else None


def fhir_datetime_or_date(value):
    if not value:
        return None
    value = value.strip()
    if FHIR_DATETIME.match(value):
        return value
    return fhir_date(value)


def with_meta(resource, meta):
    if not meta:
        return resource
    out = {}
    for key, value in resource.items():
        out[key] = value
        if key == "id":
            out["meta"] = meta
    return out


def drop_empty(resource):
    return {k: v for k, v in resource.items() if v not in (None, [], {}, "")}


def geolocation(lat, lon):
    if not (lat and lon):
        return None
    try:
        return {"url": "http://hl7.org/fhir/StructureDefinition/geolocation",
                "extension": [{"url": "latitude", "valueDecimal": float(lat)},
                              {"url": "longitude", "valueDecimal": float(lon)}]}
    except ValueError:
        return None


def administrative_code(parts):
    ext = [{"url": key, "valueCode": value} for key, value in parts if value]
    if not ext:
        return None
    return {"url": f"{KEMKES_SD}/administrativeCode", "extension": ext}


class Flags:
    def __init__(self):
        self.items = []

    def add(self, flag):
        if flag not in self.items:
            self.items.append(flag)


def create_patient(fid, sample_id, row, org, profile, flags):
    v = lambda f: get_clinical_value(row, f)  
    org_id = org.get("org_id") if org else None

    name = drop_empty({"use": "official", "family": v("family_name"),
                       "given": [v("given_name")] if v("given_name") else None})
    if len(name) == 1:
        name = None
        flags.add("patient_name_missing")

    gender_raw = (v("gender") or "").lower()
    gender = ("male" if gender_raw in ("laki-laki", "pria", "male", "m") else
              "female" if gender_raw in ("perempuan", "wanita", "female", "f") else
              "unknown" if gender_raw else None)
    if gender_raw and gender == "unknown":
        flags.add("gender_unrecognised")

    birth_date = fhir_date(v("birth_date"))
    if v("birth_date") and not birth_date:
        flags.add("birth_date_invalid_omitted")
    elif not v("birth_date"):
        flags.add("birth_date_missing")

    identifiers = []
    extensions = []
    address_ext = []
    if profile.include_nik:
        nik = v("nik")
        if nik and NIK.match(nik):
            identifiers.append({"use": "official", "system": KEMKES_NIK, "value": nik})
        elif nik:
            flags.add("nik_not_16_digits_omitted")
        else:
            flags.add("nik_missing")
        if org_id:
            identifiers.append({
                "use": "usual",
                "type": {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0203",
                                     "code": "MR", "display": "Medical record number"}]},
                "system": f"{KEMKES_SYS_IDS}/mr/{org_id}", "value": fid})
    if profile.kemkes_extensions:
        admin = administrative_code([("province", v("province_code")), ("city", v("city_code")),
                                     ("district", v("district_code")), ("village", v("village_code"))])
        if admin:
            extensions.append(admin)
            address_ext.append(admin)
        if v("citizenship_status"):
            extensions.append({"url": f"{KEMKES_SD}/citizenshipStatus",
                               "valueCode": v("citizenship_status")})
    geo = geolocation(v("latitude"), v("longitude"))
    if geo:
        address_ext.append(geo)

    country = v("country") or ("ID" if profile.kemkes_extensions else None)
    address = drop_empty({"use": "home", "type": "physical", "text": v("address"),
                          "city": v("city"), "state": v("state"), "country": country,
                          "extension": address_ext})
    has_address = any(k in address for k in ("text", "city", "state", "extension"))

    patient = {
        "resourceType": "Patient",
        "id": f"{fid}-patient",
        "active": True,
        "name": [name] if name else None,
        "gender": gender,
        "birthDate": birth_date,
        "identifier": identifiers,
        "extension": extensions,
        "address": [address] if has_address else None,
    }
    return with_meta(drop_empty(patient), profile.meta("patient"))


def create_specimen(fid, row, org, practitioner, profile, flags):
    v = lambda f: get_clinical_value(row, f)  
    org_id = org.get("org_id") if org else None
    collected = fhir_datetime_or_date(v("specimen_collected_date"))
    received = fhir_datetime_or_date(v("specimen_received_date"))
    if not collected:
        flags.add("specimen_collected_date_missing")
    if not received:
        flags.add("specimen_received_date_missing")

    collection = {}
    if collected:
        collection["collectedDateTime"] = collected
    if practitioner:
        collection["collector"] = {"reference": f"Practitioner/{practitioner['practitioner_id']}",
                                   "display": practitioner.get("name") or None}
        collection["collector"] = drop_empty(collection["collector"])
    if v("specimen_collection_method_code"):
        collection["method"] = drop_empty({
            "coding": [drop_empty({"system": "http://terminology.hl7.org/CodeSystem/v2-0488",
                                   "code": v("specimen_collection_method_code"),
                                   "display": v("specimen_collection_method_display")})],
            "text": v("specimen_collection_method_text")})
    try:
        quantity = float(v("specimen_quantity_value")) if v("specimen_quantity_value") else None
    except ValueError:
        quantity = None
        flags.add("specimen_quantity_invalid_omitted")
    if quantity is not None:
        unit = v("specimen_quantity_unit")
        collection["quantity"] = drop_empty({"value": quantity, "unit": unit,
                                             "system": "http://unitsofmeasure.org" if unit else None,
                                             "code": unit})

    specimen_type = None
    if v("specimen_type_code"):
        display = v("specimen_type_display")
        specimen_type = drop_empty({
            "coding": [drop_empty({"system": "http://snomed.info/sct", "code": v("specimen_type_code"),
                                   "display": display})],
            "text": f"{display} for TB testing" if display else None})
    else:
        flags.add("specimen_type_missing")

    identifier = profile.identifier("specimen", f"SPEC-{fid}", org_id)
    specimen = {
        "resourceType": "Specimen",
        "id": f"{fid}-specimen",
        "identifier": [identifier] if identifier else None,
        "status": "available",
        "subject": {"reference": f"Patient/{fid}-patient"},
        "receivedTime": received,
        "collection": collection,
        "type": specimen_type,
        "note": [{"text": f"Collected sputum sample ({fid}) for TB genetic testing purposes"}],
    }
    return with_meta(drop_empty(specimen), profile.meta("specimen"))


def create_encounter(fid, row, org, practitioner, profile):
    v = lambda f: get_clinical_value(row, f)  
    org_id = org.get("org_id") if org else None
    start = fhir_datetime_or_date(v("specimen_collected_date"))
    identifier = profile.identifier("encounter", fid, org_id)
    encounter = {
        "resourceType": "Encounter",
        "id": f"{fid}-encounter",
        "identifier": [identifier] if identifier else None,
        "status": "finished",
        "class": {"system": "http://terminology.hl7.org/CodeSystem/v3-ActCode", "code": "AMB",
                  "display": "ambulatory"},
        "subject": {"reference": f"Patient/{fid}-patient"},
        "period": {"start": start} if start else None,
        "serviceProvider": {"reference": f"Organization/{org_id}"} if org_id else None,
    }
    if practitioner:
        encounter["participant"] = [{
            "type": [{"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v3-ParticipationType",
                                  "code": "ATND", "display": "attender"}]}],
            "individual": drop_empty({"reference": f"Practitioner/{practitioner['practitioner_id']}",
                                      "display": practitioner.get("name")})}]
    if v("location_id"):
        encounter["location"] = [{"location": {"reference": f"Location/{v('location_id')}"}}]
    if profile.kemkes_extensions and start:
        encounter["statusHistory"] = [{"status": "finished", "period": {"start": start}}]
    return with_meta(drop_empty(encounter), profile.meta("encounter"))


def create_service_request(fid, org, practitioner, profile, timestamp):
    org_id = org.get("org_id") if org else None
    identifier = profile.identifier("service_request", f"SR-{fid}", org_id)
    request = {
        "resourceType": "ServiceRequest",
        "id": f"{fid}-service-request",
        "identifier": [identifier] if identifier else None,
        "status": "active",
        "intent": "original-order",
        "priority": "routine",
        "category": [{"coding": [{"system": "http://snomed.info/sct", "code": "108252007",
                                  "display": "Laboratory procedure"}]}],
        "code": {"coding": [{"system": "http://loinc.org", "code": "69548-6",
                             "display": "Genetic variant assessment"}],
                 "text": "TB Genetic Variant Assessment"},
        "subject": {"reference": f"Patient/{fid}-patient"},
        "encounter": {"reference": f"Encounter/{fid}-encounter", "display": "TB Testing Encounter"},
        "occurrenceDateTime": timestamp,
    }
    if practitioner:
        request["requester"] = drop_empty({"reference": f"Practitioner/{practitioner['practitioner_id']}",
                                           "display": practitioner.get("name")})
        if practitioner.get("role_id"):
            request["performer"] = [drop_empty({"reference": f"PractitionerRole/{practitioner['role_id']}",
                                                "display": practitioner.get("role_display")})]
    return with_meta(drop_empty(request), profile.meta("service_request"))


def create_organization(org, profile):
    g = lambda f: org.get(f) or None  
    telecom = [t for t in (
        {"system": "phone", "value": g("phone"), "use": "work"} if g("phone") else None,
        {"system": "email", "value": g("email"), "use": "work"} if g("email") else None) if t]
    address_ext = []
    if profile.kemkes_extensions:
        admin = administrative_code([("province", g("province_code")), ("city", g("city_code")),
                                     ("district", g("district_code"))])
        if admin:
            address_ext.append(admin)
    geo = geolocation(g("latitude"), g("longitude"))
    if geo:
        address_ext.append(geo)
    country = g("country") or ("ID" if profile.kemkes_extensions else None)
    address = drop_empty({"use": "work", "type": "physical",
                          "line": [g("address_line")] if g("address_line") else None,
                          "city": g("city"), "state": g("state"), "country": country,
                          "extension": address_ext})

    org_type = None
    if profile.kemkes_extensions and g("type_code"):
        org_type = drop_empty({"coding": [drop_empty({
            "system": "http://terminology.kemkes.go.id/CodeSystem/organization-type",
            "code": g("type_code"), "display": g("type_display")})], "text": g("type_text")})
    elif g("type_text") or g("type_display"):
        org_type = {"text": g("type_text") or g("type_display")}

    identifier = profile.identifier("organization", org["org_id"])
    resource = {
        "resourceType": "Organization",
        "id": org["org_id"],
        "identifier": [identifier] if identifier else None,
        "active": True,
        "type": [org_type] if org_type else None,
        "name": g("name"),
        "alias": [g("alias")] if g("alias") else None,
        "telecom": telecom,
        "address": [address] if len(address) > 2 else None,
    }
    return with_meta(drop_empty(resource), profile.meta("organization"))


def create_practitioner(practitioner, profile):
    g = lambda f: practitioner.get(f) or None  
    resource = {
        "resourceType": "Practitioner",
        "id": practitioner["practitioner_id"],
        "active": True,
        "name": [{"use": "official", "text": g("name")}] if g("name") else None,
        "telecom": [{"system": "phone", "value": g("phone"), "use": "work"}] if g("phone") else None,
        "gender": g("gender"),
        "birthDate": fhir_date(g("birth_date")),
    }
    if profile.include_nik and g("nik"):
        resource["identifier"] = [{"use": "official", "system": KEMKES_NIK, "value": g("nik")}]
    if profile.kemkes_extensions and g("str_kki_number"):
        qualification = {
            "identifier": [{"system": "https://fhir.kemkes.go.id/id/str-kki-number",
                            "value": g("str_kki_number")}],
            "code": {"coding": [{"system": "https://terminology.kemkes.go.id/v1-0302",
                                 "code": "STR-KKI", "display": "Surat Tanda Registrasi Dokter"}],
                     "text": "Surat Tanda Registrasi Dokter"},
        }
        if g("qualification_period_start"):
            qualification["period"] = {"start": g("qualification_period_start")}
        resource["qualification"] = [qualification]
    return with_meta(drop_empty(resource), profile.meta("practitioner"))


def create_practitioner_role(practitioner, org, profile):
    g = lambda f: practitioner.get(f) or None 
    if not g("role_id"):
        return None
    resource = {
        "resourceType": "PractitionerRole",
        "id": g("role_id"),
        "active": True,
        "practitioner": drop_empty({"reference": f"Practitioner/{practitioner['practitioner_id']}",
                                    "display": g("name")}),
        "organization": drop_empty({"reference": f"Organization/{org['org_id']}",
                                    "display": org.get("name")}) if org else None,
        "code": [drop_empty({"coding": [drop_empty({"system": "http://snomed.info/sct",
                                                    "code": g("role_code"),
                                                    "display": g("role_display")})],
                             "text": g("role_text")})] if g("role_code") else None,
        "telecom": [{"system": "phone", "value": g("phone"), "use": "work"}] if g("phone") else None,
    }
    return with_meta(drop_empty(resource), profile.meta("practitioner_role"))


def collection(resources, *id_parts):
    return {"resourceType": "Bundle", "id": bundle_id("clinical", *id_parts), "type": "collection",
            "timestamp": now_iso(),
            "entry": [{"fullUrl": full_url(r["resourceType"], r["id"]), "resource": r} for r in resources]}


def read_sample_ids(path):
    with open(path, encoding="utf-8") as handle:
        return [ln.strip() for ln in handle if ln.strip() and not ln.startswith("#")]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--patient_metadata', required=True)
    p.add_argument('--organization_metadata')
    p.add_argument('--practitioner_metadata')
    p.add_argument('--sample_ids', required=True, help='file with one sample id per line')
    p.add_argument('--fhir_profile', default='generic', choices=['generic', 'satusehat'])
    p.add_argument('--identifier_system_base', default=None)
    p.add_argument('--outdir', default='.')
    args = p.parse_args(argv)

    profile = get_profile(args.fhir_profile, args.identifier_system_base)
    try:
        patients = load_clinical_metadata(args.patient_metadata)
        org = load_organization_metadata(args.organization_metadata) if args.organization_metadata else None
        practitioner = (load_practitioner_metadata(args.practitioner_metadata)
                        if args.practitioner_metadata else None)
    except MetadataError as exc:
        sys.exit(f"ERROR: {exc}")
    sample_ids = read_sample_ids(args.sample_ids)
    os.makedirs(args.outdir, exist_ok=True)
    timestamp = now_iso()

    shared = []
    if org:
        shared.append(create_organization(org, profile))
    if practitioner:
        shared.append(create_practitioner(practitioner, profile))
        role = create_practitioner_role(practitioner, org, profile)
        if role:
            shared.append(role)
    with open(os.path.join(args.outdir, "shared.clinical.json"), "w", encoding="utf-8") as out:
        json.dump(collection(shared, "shared"), out, indent=2, ensure_ascii=False)

    summary = []
    external = {}
    used_rows = set()
    for sample_id in sample_ids:
        fid = fhir_id(sample_id)
        row = lookup_sample(patients, sample_id)
        if row is None:
            summary.append((sample_id, fid, "missing_clinical", ""))
            continue
        used_rows.add(row["sample_id"])
        flags = Flags()
        resources = [create_patient(fid, sample_id, row, org, profile, flags),
                     create_specimen(fid, row, org, practitioner, profile, flags)]
        encounter_id = get_clinical_value(row, "encounter_id")
        if encounter_id:
            external.setdefault(fid, {})["Encounter"] = encounter_id
        else:
            resources.append(create_encounter(fid, row, org, practitioner, profile))
        resources.append(create_service_request(fid, org, practitioner, profile, timestamp))
        with open(os.path.join(args.outdir, f"{fid}.clinical.json"), "w", encoding="utf-8") as out:
            json.dump(collection(resources, fid), out, indent=2, ensure_ascii=False)
        summary.append((sample_id, fid, "matched", ",".join(flags.items)))

    for sample_id in sorted(set(patients) - used_rows):
        summary.append((sample_id, fhir_id(sample_id), "orphan_row", ""))

    with open(os.path.join(args.outdir, "external_refs.json"), "w", encoding="utf-8") as out:
        json.dump(external, out, indent=2, sort_keys=True)
    with open(os.path.join(args.outdir, "clinical_summary.tsv"), "w", encoding="utf-8") as out:
        out.write("sample_id\tfhir_id\tstatus\tflags\n")
        for row in summary:
            out.write("\t".join(row) + "\n")

    matched = sum(1 for r in summary if r[2] == "matched")
    missing = [r[0] for r in summary if r[2] == "missing_clinical"]
    orphans = [r[0] for r in summary if r[2] == "orphan_row"]
    print(f"clinical resources for {matched}/{len(sample_ids)} sample(s) (profile {profile.name})")
    if missing:
        print(f"WARNING: no clinical row for: {', '.join(missing)}")
    if orphans:
        print(f"WARNING: clinical rows with no genomic sample: {', '.join(orphans)}")


if __name__ == "__main__":
    main()
