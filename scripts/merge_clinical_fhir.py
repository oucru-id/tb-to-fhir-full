#!/usr/bin/env python3

import json
import argparse
import uuid
import sys
import os  
from datetime import datetime, timezone
from clinical_metadata_parser import (
    load_clinical_metadata, find_matching_sample, get_clinical_value,
    load_organization_metadata, load_practitioner_metadata
)
import base64
import re

def debug_print(message):
    print(f"DEBUG: {message}", file=sys.stderr)

def create_patient_resource(sample_id, clinical_data=None, org_data=None):
    if not clinical_data:
        raise ValueError(f"Clinical data is required for sample {sample_id}")

    org_data = org_data or {}
    org_id = org_data.get('org_id')

    family_name = get_clinical_value(clinical_data, 'family_name')
    given_name = get_clinical_value(clinical_data, 'given_name')
    gender = get_clinical_value(clinical_data, 'gender', 'unknown').lower()
    birth_date = get_clinical_value(clinical_data, 'birth_date')
    nik = get_clinical_value(clinical_data, 'nik')
    province_code = get_clinical_value(clinical_data, 'province_code')
    city_code = get_clinical_value(clinical_data, 'city_code')
    district_code = get_clinical_value(clinical_data, 'district_code')
    village_code = get_clinical_value(clinical_data, 'village_code')
    citizenship_status = get_clinical_value(clinical_data, 'citizenship_status')
    lat = get_clinical_value(clinical_data, 'latitude', None)
    lon = get_clinical_value(clinical_data, 'longitude', None)

    if gender in ['laki-laki', 'pria', 'male', 'm']:
        gender = "male"
    elif gender in ['perempuan', 'wanita', 'female', 'f']:
        gender = "female"
    else:
        gender = "unknown"

    geo_extensions = []
    if lat and lon:
        try:
            geo_extensions = [{
                "url": "http://hl7.org/fhir/StructureDefinition/geolocation",
                "extension": [
                    {"url": "latitude",  "valueDecimal": float(lat)},
                    {"url": "longitude", "valueDecimal": float(lon)}
                ]
            }]
        except ValueError:
            pass

    return {
        "resourceType": "Patient",
        "id": f"{sample_id}-patient",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Patient"]
        },
        "active": True,
        "name": [
            {
                "use": "official",
                "family": family_name,
                "given": [given_name]
            }
        ],
        "gender": gender,
        "birthDate": birth_date,
        "identifier": [
            {
                "use": "official",
                "system": "https://fhir.kemkes.go.id/id/nik",
                "value": nik
            },
            {
                "use": "usual",
                "type": {
                    "coding": [
                        {
                            "system": "http://terminology.hl7.org/CodeSystem/v2-0203",
                            "code": "MR",
                            "display": "Medical record number"
                        }
                    ]
                },
                "system": f"http://sys-ids.kemkes.go.id/mr/{org_id}",
                "value": sample_id
            }
        ],
        "extension": [
            {
                "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode",
                "extension": [
                    {"url": "province", "valueCode": province_code},
                    {"url": "city",     "valueCode": city_code},
                    {"url": "district", "valueCode": district_code},
                    {"url": "village",  "valueCode": village_code}
                ]
            },
            {
                "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/citizenshipStatus",
                "valueCode": citizenship_status
            }
        ],
        "address": [
            {
                "use": "home",
                "type": "physical",
                "text": get_clinical_value(clinical_data, 'address'),
                "city": get_clinical_value(clinical_data, 'city'),
                "state": get_clinical_value(clinical_data, 'state'),
                "country": "ID",
                "extension": [
                    {
                        "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode",
                        "extension": [
                            {"url": "province", "valueCode": province_code},
                            {"url": "city",     "valueCode": city_code},
                            {"url": "district", "valueCode": district_code},
                            {"url": "village",  "valueCode": village_code}
                        ]
                    },
                    *geo_extensions
                ]
            }
        ]
    }

def create_specimen_resource(sample_id, clinical_data=None, practitioner_data=None, org_data=None):
    org_data = org_data or {}
    practitioner_data = practitioner_data or {}

    org_id = org_data.get('org_id',)
    practitioner_id = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    practitioner_name = practitioner_data.get('name', 'Unknown Practitioner')

    if clinical_data:
        given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
        spec_type_code    = get_clinical_value(clinical_data, 'specimen_type_code', '119334006')
        spec_type_display = get_clinical_value(clinical_data, 'specimen_type_display', 'Sputum specimen')
        method_code       = get_clinical_value(clinical_data, 'specimen_collection_method_code', 'KOFFP')
        method_display    = get_clinical_value(clinical_data, 'specimen_collection_method_display', 'Plate, Cough')
        method_text       = get_clinical_value(clinical_data, 'specimen_collection_method_text', 'Sputum collection')
        qty_value         = get_clinical_value(clinical_data, 'specimen_quantity_value', '2')
        qty_unit          = get_clinical_value(clinical_data, 'specimen_quantity_unit', 'mL')
        collected_date    = get_clinical_value(clinical_data, 'specimen_collected_date', None)
        received_date     = get_clinical_value(clinical_data, 'specimen_received_date', None)
    else:
        patient_display   = f"Patient {sample_id}"
        spec_type_code    = '119334006'
        spec_type_display = 'Sputum specimen'
        method_code       = 'KOFFP'
        method_display    = 'Plate, Cough'
        method_text       = 'Sputum collection'
        qty_value         = '2'
        qty_unit          = 'mL'
        collected_date    = None
        received_date     = None

    now = datetime.now(timezone.utc).isoformat()
    collected_dt = collected_date if collected_date and collected_date != 'Unknown' else now
    received_dt  = received_date  if received_date  and received_date  != 'Unknown' else now

    try:
        qty_float = float(qty_value)
    except (ValueError, TypeError):
        qty_float = 2.0

    return {
        "resourceType": "Specimen",
        "id": f"{sample_id}-specimen",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Specimen"]
        },
        "identifier": [
            {
                "system": f"http://sys-ids.kemkes.go.id/specimen/{org_id}",
                "value": f"SPEC-{sample_id}"
            }
        ],
        "status": "available",
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "receivedTime": received_dt,
        "collection": {
            "collectedDateTime": collected_dt,
            "collector": {
                "reference": f"Practitioner/{practitioner_id}",
                "display": practitioner_name
            },
            "method": {
                "coding": [
                    {
                        "system": "http://terminology.hl7.org/CodeSystem/v2-0488",
                        "code": method_code,
                        "display": method_display
                    }
                ],
                "text": method_text
            },
            "quantity": {
                "value": qty_float,
                "unit": qty_unit,
                "system": "http://unitsofmeasure.org",
                "code": qty_unit
            }
        },
        "type": {
            "coding": [
                {
                    "system": "http://snomed.info/sct",
                    "code": spec_type_code,
                    "display": spec_type_display
                }
            ],
            "text": f"{spec_type_display} for TB testing"
        },
        "note": [
            {
                "text": f"Collected sputum sample from {patient_display} ({sample_id}) for TB genetic testing purposes"
            }
        ]
    }

def create_organization_resource(org_data=None):
    org_data = org_data or {}
    org_id       = org_data.get('org_id', 'unknown-org')
    name         = org_data.get('name', 'Unknown Organization')
    alias        = org_data.get('alias', '')
    type_code    = org_data.get('type_code', '')
    type_display = org_data.get('type_display', '')
    type_text    = org_data.get('type_text', '')
    phone        = org_data.get('phone', '')
    email        = org_data.get('email', '')
    address_line = org_data.get('address_line', '')
    city         = org_data.get('city', '')
    state        = org_data.get('state', '')
    country      = org_data.get('country', 'ID')
    province_code = org_data.get('province_code', '')
    city_code     = org_data.get('city_code', '')
    district_code = org_data.get('district_code', '')
    lat           = org_data.get('latitude', None)
    lon           = org_data.get('longitude', None)

    telecom = []
    if phone:
        telecom.append({"system": "phone", "value": phone, "use": "work"})
    if email:
        telecom.append({"system": "email", "value": email, "use": "work"})

    addr_extensions = []
    if province_code or city_code or district_code:
        code_ext = {"url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode", "extension": []}
        if province_code:
            code_ext["extension"].append({"url": "province", "valueCode": province_code})
        if city_code:
            code_ext["extension"].append({"url": "city",     "valueCode": city_code})
        if district_code:
            code_ext["extension"].append({"url": "district", "valueCode": district_code})
        addr_extensions.append(code_ext)
    if lat and lon:
        try:
            addr_extensions.append({
                "url": "http://hl7.org/fhir/StructureDefinition/geolocation",
                "extension": [
                    {"url": "latitude",  "valueDecimal": float(lat)},
                    {"url": "longitude", "valueDecimal": float(lon)}
                ]
            })
        except ValueError:
            pass

    resource = {
        "resourceType": "Organization",
        "id": org_id,
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Organization"]
        },
        "identifier": [
            {
                "use": "official",
                "system": "http://sys-ids.kemkes.go.id/organization",
                "value": org_id
            }
        ],
        "active": True,
        "type": [
            {
                "coding": [
                    {
                        "system": "http://terminology.kemkes.go.id/CodeSystem/organization-type",
                        "code": type_code,
                        "display": type_display
                    }
                ],
                "text": type_text
            }
        ],
        "name": name,
        "telecom": telecom,
        "address": [
            {
                "use": "work",
                "type": "physical",
                "line": [address_line] if address_line else [],
                "city": city,
                "state": state,
                "country": country,
                "extension": addr_extensions
            }
        ]
    }
    if alias:
        resource["alias"] = [alias]
    return resource

def create_practitioner_resource(practitioner_data=None):
    practitioner_data = practitioner_data or {}
    pid       = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    nik       = practitioner_data.get('nik', '')
    name      = practitioner_data.get('name', 'Unknown Practitioner')
    phone     = practitioner_data.get('phone', '')
    gender    = practitioner_data.get('gender', 'unknown')
    birth_date = practitioner_data.get('birth_date', '')
    str_kki   = practitioner_data.get('str_kki_number', '')
    qual_start = practitioner_data.get('qualification_period_start', '')

    telecom = []
    if phone:
        telecom.append({"system": "phone", "value": phone, "use": "work"})

    qualification = []
    if str_kki:
        qual = {
            "code": {
                "coding": [{
                    "system": "https://terminology.kemkes.go.id/v1-0302",
                    "code": "STR-KKI",
                    "display": "Surat Tanda Registrasi Dokter"
                }],
                "text": "Surat Tanda Registrasi Dokter"
            }
        }
        if str_kki:
            qual["identifier"] = [{"system": "https://fhir.kemkes.go.id/id/str-kki-number", "value": str_kki}]
        if qual_start:
            qual["period"] = {"start": qual_start}
        qualification.append(qual)

    resource = {
        "resourceType": "Practitioner",
        "id": pid,
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Practitioner"]
        },
        "active": True,
        "name": [{"use": "official", "text": name}],
        "telecom": telecom,
        "gender": gender
    }
    if nik:
        resource["identifier"] = [{
            "use": "official",
            "system": "https://fhir.kemkes.go.id/id/nik",
            "value": nik
        }]
    if birth_date:
        resource["birthDate"] = birth_date
    if qualification:
        resource["qualification"] = qualification
    return resource

def create_practitioner_role_resource(practitioner_data=None, org_data=None):
    practitioner_data = practitioner_data or {}
    org_data = org_data or {}
    pid          = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    pname        = practitioner_data.get('name', 'Unknown Practitioner')
    phone        = practitioner_data.get('phone', '')
    role_id      = practitioner_data.get('role_id', 'unknown-role')
    role_code    = practitioner_data.get('role_code', '')
    role_display = practitioner_data.get('role_display', '')
    role_text    = practitioner_data.get('role_text', '')
    org_id       = org_data.get('org_id', 'unknown-org')
    org_name     = org_data.get('name', 'Unknown Organization')

    telecom = []
    if phone:
        telecom.append({"system": "phone", "value": phone, "use": "work"})

    resource = {
        "resourceType": "PractitionerRole",
        "id": role_id,
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/PractitionerRole"]
        },
        "active": True,
        "practitioner": {
            "reference": f"Practitioner/{pid}",
            "display": pname
        },
        "organization": {
            "reference": f"Organization/{org_id}",
            "display": org_name
        },
        "code": [
            {
                "coding": [
                    {
                        "system": "http://snomed.info/sct",
                        "code": role_code,
                        "display": role_display
                    }
                ],
                "text": role_text
            }
        ],
        "telecom": telecom
    }
    return resource

def classify_drug_resistance(observations):
    resistant_drugs_groups = set()
    detected_drugs = set()
    detected_genes = set()
    
    for obs in observations:
        codes = obs.get('code', {}).get('coding', [])
        is_panel = any(c.get('code') == '89486-5' for c in codes)
        
        if is_panel:
            components = obs.get('component', [])
            for comp in components:
                value_coding = comp.get('valueCodeableConcept', {}).get('coding', [])
                is_resistant = any(vc.get('code') == 'LA6676-6' for vc in value_coding)
                
                if is_resistant:
                    drug_display = comp.get('code', {}).get('coding', [{}])[0].get('display', '').lower()
                    
                    detected_drugs.add(drug_display.replace('[susceptibility]', '').strip())
                    
                    if 'rifampicin' in drug_display or 'rifampin' in drug_display:
                        resistant_drugs_groups.add('rifampicin')
                    elif 'isoniazid' in drug_display:
                        resistant_drugs_groups.add('isoniazid')
                    elif 'ethambutol' in drug_display:
                        resistant_drugs_groups.add('ethambutol')
                    elif 'pyrazinamide' in drug_display:
                        resistant_drugs_groups.add('pyrazinamide')
                    elif 'streptomycin' in drug_display:
                        resistant_drugs_groups.add('streptomycin')
                    elif 'ethionamide' in drug_display:
                        resistant_drugs_groups.add('ethionamide')
                    elif any(fq in drug_display for fq in ['levofloxacin', 'moxifloxacin', 'ofloxacin', 'ciprofloxacin']):
                        resistant_drugs_groups.add('fluoroquinolone')
                    elif any(sli in drug_display for sli in ['amikacin', 'kanamycin', 'capreomycin']):
                        resistant_drugs_groups.add('second_line_injectable')
                    elif any(ga in drug_display for ga in ['bedaquiline', 'linezolid']):
                        resistant_drugs_groups.add('group_a')
                    else:
                        resistant_drugs_groups.add(drug_display)

    for obs in observations:
        codes = obs.get('code', {}).get('coding', [])
        is_variant = any(c.get('code') == '69548-6' for c in codes)
        
        if is_variant:
            components = obs.get('component', [])
            
            significance = ""
            current_gene = ""
            
            for component in components:
                code_display = component.get('code', {}).get('coding', [{}])[0].get('display', '').lower()
                
                if 'gene studied' in code_display:
                    current_gene = component.get('valueCodeableConcept', {}).get('text', '')
                
                if 'genetic variation clinical significance' in code_display:
                    significance = component.get('valueCodeableConcept', {}).get('text', '')

            if current_gene and significance and "Assoc w R" in significance:
                detected_genes.add(current_gene)
    
    has_rif = 'rifampicin' in resistant_drugs_groups
    has_inh = 'isoniazid' in resistant_drugs_groups
    has_fq = 'fluoroquinolone' in resistant_drugs_groups
    has_group_a = 'group_a' in resistant_drugs_groups
    
    is_mdr = has_rif and has_inh
    is_rr = has_rif
    
    classification = "Sensitive"
    description = "No resistance mutations detected"

    if not resistant_drugs_groups:
        classification = "Sensitive"
        description = "No resistance mutations detected"
    elif (is_mdr or is_rr) and has_fq and has_group_a:
        classification = "XDR-TB"
        description = "Extensively drug-resistant tuberculosis (MDR/RR + FQ + Group A)"
    elif (is_mdr or is_rr) and has_fq:
        classification = "Pre-XDR-TB"
        description = "Pre-extensively drug-resistant tuberculosis (MDR/RR + FQ)"
    elif has_rif and has_inh:
        classification = "MDR-TB"
        description = "Multidrug-resistant tuberculosis"
    elif has_rif and not has_inh:
        classification = "RR-TB"
        description = "Rifampicin-resistant tuberculosis"
    elif has_inh and not has_rif:
        classification = "HR-TB"
        description = "Isoniazid-resistant tuberculosis"
    elif len(resistant_drugs_groups) == 1:
        if 'streptomycin' in resistant_drugs_groups:
            classification = "Streptomycin-resistant TB"
            description = "Streptomycin mono-resistant tuberculosis"
        elif 'ethionamide' in resistant_drugs_groups:
            classification = "Ethionamide-resistant TB"
            description = "Ethionamide mono-resistant tuberculosis"
        elif 'pyrazinamide' in resistant_drugs_groups:
            classification = "Pyrazinamide-resistant TB"
            description = "Pyrazinamide mono-resistant tuberculosis"
        elif 'ethambutol' in resistant_drugs_groups:
            classification = "Ethambutol-resistant TB"
            description = "Ethambutol mono-resistant tuberculosis"
        elif 'fluoroquinolone' in resistant_drugs_groups:
            has_cipro = any('ciprofloxacin' in d for d in detected_drugs)
            has_other_fq = any(fq in d for d in detected_drugs for fq in ['levofloxacin', 'moxifloxacin', 'ofloxacin'])
            
            if has_cipro and not has_other_fq:
                classification = "Ciprofloxacin-resistant TB"
                description = "Ciprofloxacin mono-resistant tuberculosis"
            else:
                classification = "Drug-resistant"
                description = f"Resistance to: {', '.join(sorted(resistant_drugs_groups))}"
        else:
            classification = "Drug-resistant"
            description = f"Resistance to: {', '.join(sorted(resistant_drugs_groups))}"
    else:
        classification = "Drug-resistant"
        description = f"Resistance to: {', '.join(sorted(resistant_drugs_groups))}"

    return classification, description, sorted(list(detected_genes)), sorted(list(detected_drugs))

def extract_lineage_info(observations):
    for obs in observations:
        codes = obs.get('code', {}).get('coding', [])
        for code in codes:
            if code.get('code') == '614-8':
                lineage_text = obs.get('valueCodeableConcept', {}).get('text', '')
                if lineage_text and lineage_text != 'unknown':
                    return lineage_text

        components = obs.get('component', [])
        for component in components:
            code_display = component.get('code', {}).get('coding', [{}])[0].get('display', '').lower()
            if 'mycobacterial strain' in code_display or 'lineage' in code_display:
                lineage_text = component.get('valueCodeableConcept', {}).get('text', '')
                if lineage_text and lineage_text != 'unknown':
                    return lineage_text
    return None

def get_resistance_conclusion_coding(resistance_class):
    coding_map = {
        "Sensitive": {
           "system": "https://terminology.kemkes.go.id/CodeSystem/episodeofcare-type",
           "code": "TB-SO",
           "display": "Tuberkulosis Sensitif Obat"
        },
        "RR-TB": {
            "system": "http://snomed.info/sct",
            "code": "415345001",
            "display": "Rifampicin resistant tuberculosis"
        },
        "HR-TB": {
            "system": "http://snomed.info/sct", 
            "code": "414546009",
            "display": "Isoniazid resistant tuberculosis"
        },
        "MDR-TB": {
            "system": "http://snomed.info/sct",
            "code": "423092005", 
            "display": "Multidrug resistant tuberculosis"
        },
        "Pre-XDR-TB": {
            "system": "http://terminology.kemkes.go.id/CodeSystem/clinical-term",
            "code": "OV000435",
            "display": "Pre-XDR"
        },
        "XDR-TB": {
            "system": "http://snomed.info/sct",
            "code": "710106005",
            "display": "Extensively drug resistant tuberculosis"
        },
        "Streptomycin-resistant TB": {
            "system": "http://snomed.info/sct",
            "code": "415622003",
            "display": "Streptomycin resistant tuberculosis"
        },
        "Ethionamide-resistant TB": {
            "system": "http://snomed.info/sct",
            "code": "414149006",
            "display": "Ethionamide resistant tuberculosis"
        },
        "Pyrazinamide-resistant TB": {
            "system": "http://snomed.info/sct",
            "code": "415222009",
            "display": "Pyrazinamide resistant tuberculosis"
        },
        "Ciprofloxacin-resistant TB": {
            "system": "http://snomed.info/sct",
            "code": "413852006",
            "display": "Ciprofloxacin resistant tuberculosis"
        },
        "Ethambutol-resistant TB": {
            "system": "http://snomed.info/sct",
            "code": "414146004",
            "display": "Ethambutol resistant tuberculosis"
        },
        "Drug-resistant": {
            "system": "http://snomed.info/sct",
            "code": "413556004",
            "display": "Antibiotic resistant tuberculosis"
        }
    }
    
    return coding_map.get(resistance_class)

def create_diagnostic_report(sample_id, observations, clinical_data=None, org_data=None, practitioner_data=None):
    
    resistance_class, resistance_description, resistant_genes, resistant_drugs = classify_drug_resistance(observations)
    
    lineage_info = extract_lineage_info(observations)
    
    conclusion_parts = []
    
    if resistance_class == "Sensitive":
        conclusion_parts.append(resistance_class)
        conclusion_parts.append("No resistance-associated gene detected")
        conclusion_parts.append("No drug resistance detected")
    else:
        conclusion_parts.append(f"{resistance_class} ({resistance_description})")
        if resistant_genes:
            conclusion_parts.append(f"Detected resistance genes: {', '.join(resistant_genes)}")
        if resistant_drugs:
            conclusion_parts.append(f"Detected drug resistance: {', '.join(resistant_drugs)}")
    
    if lineage_info:
        conclusion_parts.append(f"TB {lineage_info} detected")
    
    conclusion_parts.append("Reference genome: NC_000962.3")
    
    conclusion = ". ".join(conclusion_parts)
    
    conclusion_codes = []
    
    if resistance_class:
        resistance_coding = get_resistance_conclusion_coding(resistance_class)
        
        display_text = resistance_class
        if resistance_class == "Sensitive":
            display_text = "Sensitive - No resistance detected"
            
        code_entry = {"text": display_text}
        
        if resistance_coding:
            code_entry["coding"] = [resistance_coding]
            
        conclusion_codes.append(code_entry)
    
    if lineage_info:
        conclusion_codes.append({
            "text": lineage_info
        })
    
    org_data = org_data or {}
    practitioner_data = practitioner_data or {}
    org_id   = org_data.get('org_id', 'unknown-org')
    org_name = org_data.get('name', 'Unknown Organization')
    pid      = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    pname    = practitioner_data.get('name', 'Unknown Practitioner')

    if clinical_data:
        given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
    else:
        patient_display = f"Patient {sample_id}"

    report_id = f"{sample_id}-genomic-report"
    current_time = datetime.now(timezone.utc).isoformat()
    
    html_content = f"""<div xmlns="http://www.w3.org/1999/xhtml">
<h1>TB Genomic Analysis Report</h1>
<p><strong>Patient:</strong> {patient_display}</p>
<p><strong>Sample ID:</strong> {sample_id}</p>
<p><strong>Report Date:</strong> {current_time}</p>
<p><strong>Resistance Classification:</strong> {resistance_class}</p>
<p><strong>Conclusion:</strong> {conclusion}</p>
"""
    
    if observations:
        html_content += "<h2>Detected Mutations</h2><ul>"
        for obs in observations:
            components = obs.get('component', [])
            for component in components:
                if 'valueCodeableConcept' in component:
                    text = component['valueCodeableConcept'].get('text', '')
                    if text and text != 'unknown':
                        html_content += f"<li>{text}</li>"
        html_content += "</ul>"
    
    if lineage_info:
        html_content += f"<p><strong>Mycobacterial Lineage:</strong> {lineage_info}</p>"
    
    html_content += "</div>"
    
    html_base64 = base64.b64encode(html_content.encode('utf-8')).decode('utf-8')
    
    return {
        "resourceType": "DiagnosticReport",
        "id": report_id,
        "meta": {
            "profile": ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomics-report"],
            "tag": [
                {
                    "system": "http://terminology.kemkes.go.id/sp",
                    "code": "genomics",
                    "display": "Genomics"
                }
            ]
        },
        "identifier": [
            {
                "system": f"http://sys-ids.kemkes.go.id/diagnostic-report/{org_id}",
                "value": f"TB-GEN-{sample_id}-{datetime.now().strftime('%Y%m%d')}"
            }
        ],
        "status": "final",
        "category": [
            {
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                    "code": "GE",
                    "display": "Genetics"
                }]
            }
        ],
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": "81247-9",
                "display": "Master HL7 genetic variant reporting panel"
            }],
            "text": "TB Genomic Analysis Report"
        },
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "encounter": {
            "reference": f"Encounter/{sample_id}-encounter",
            "display": "TB Testing Encounter"
        },
        "effectiveDateTime": current_time,
        "issued": current_time,
        "performer": [
            {
                "reference": f"Organization/{org_id}",
                "display": org_name
            },
            {
                "reference": f"Practitioner/{pid}",
                "display": pname
            }
        ],
        "result": [{"reference": f"Observation/{obs['id']}"} for obs in observations if obs.get('id')],
        "specimen": [{
            "reference": f"Specimen/{sample_id}-specimen",
            "display": f"Sputum specimen from {patient_display}"
        }],
        "conclusion": conclusion,
        "conclusionCode": conclusion_codes,
        "presentedForm": [
            {
                "contentType": "text/html",
                "language": "en-US", 
                "title": "TB Genomic Analysis Report",
                "data": html_base64
            }
        ]
    }

def create_service_request_resource(sample_id, clinical_data=None, practitioner_data=None, org_data=None):
    practitioner_data = practitioner_data or {}
    org_data = org_data or {}
    org_id       = org_data.get('org_id', 'unknown-org')
    pid          = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    pname        = practitioner_data.get('name', 'Unknown Practitioner')
    role_id      = practitioner_data.get('role_id', 'unknown-role')
    role_display = practitioner_data.get('role_display', 'Unknown role')

    if clinical_data:
        given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
    else:
        patient_display = f"Patient {sample_id}"

    return {
        "resourceType": "ServiceRequest",
        "id": f"{sample_id}-service-request",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/ServiceRequest"]
        },
        "identifier": [
            {
                "system": f"http://sys-ids.kemkes.go.id/servicerequest/{org_id}",
                "value": f"SR-{sample_id}"
            }
        ],
        "status": "active",
        "intent": "original-order",
        "priority": "routine",
        "category": [
            {
                "coding": [
                    {
                        "system": "http://snomed.info/sct",
                        "code": "108252007",
                        "display": "Laboratory procedure"
                    }
                ]
            }
        ],
        "code": {
            "coding": [
                {
                    "system": "http://loinc.org",
                    "code": "69548-6",
                    "display": "Genetic variant assessment"
                }
            ],
            "text": "TB Genetic Variant Assessment"
        },
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "encounter": {
            "reference": f"Encounter/{sample_id}-encounter",
            "display": "TB Testing Encounter"
        },
        "occurrenceDateTime": datetime.now(timezone.utc).isoformat(),
        "requester": {
            "reference": f"Practitioner/{pid}",
            "display": pname
        },
        "performer": [
            {
                "reference": f"PractitionerRole/{role_id}",
                "display": role_display
            }
        ]
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Path to input FHIR bundle')
    parser.add_argument('--output', required=True, help='Path to output merged FHIR bundle')
    parser.add_argument('--patient_metadata',      help='Path to patient_clinical_metadata CSV/Excel file')
    parser.add_argument('--organization_metadata', help='Path to organization_metadata CSV/Excel file')
    parser.add_argument('--practitioner_metadata', help='Path to practitioner_metadata CSV/Excel file')
    args = parser.parse_args()

    debug_print(f"Input FHIR file: {args.input}")
    debug_print(f"Patient metadata: {args.patient_metadata}")
    debug_print(f"Organization metadata: {args.organization_metadata}")
    debug_print(f"Practitioner metadata: {args.practitioner_metadata}")
    debug_print(f"Output file: {args.output}")

    clinical_data = {}
    if args.patient_metadata and os.path.exists(args.patient_metadata):
        clinical_data = load_clinical_metadata(args.patient_metadata)
    else:
        debug_print(f"Patient metadata file not found or not provided: {args.patient_metadata}")

    org_data = {}
    if args.organization_metadata and os.path.exists(args.organization_metadata):
        org_data = load_organization_metadata(args.organization_metadata)
    else:
        debug_print(f"Organization metadata file not found or not provided: {args.organization_metadata}")

    practitioner_data = {}
    if args.practitioner_metadata and os.path.exists(args.practitioner_metadata):
        practitioner_data = load_practitioner_metadata(args.practitioner_metadata)
    else:
        debug_print(f"Practitioner metadata file not found or not provided: {args.practitioner_metadata}")

    try:
        with open(args.input, 'r') as f:
            fhir_bundle = json.load(f)

        sample_ids = set()
        all_observations = []
        
        for entry in fhir_bundle.get('entry', []):
            resource = entry.get('resource', {})
            if resource.get('resourceType') == 'Observation':
                all_observations.append(resource)
                subject_ref = resource.get('subject', {}).get('reference', '')
                
                if subject_ref.startswith('Patient/'):
                    sample_id = subject_ref.replace('Patient/', '').replace('-patient', '')
                    sample_ids.add(sample_id)

        if not sample_ids or all(sid.startswith('NC-') for sid in sample_ids):
            filename = os.path.basename(args.input)
            filename_sample_id = filename.replace('.fhir.json', '').replace('_ont', '').replace('_illumina', '')
            sample_ids.add(filename_sample_id)

        matched_samples = {}
        for sample_id in sample_ids:
            sample_clinical_data = find_matching_sample(sample_id, clinical_data)
            if sample_clinical_data:
                matched_samples[sample_id] = sample_clinical_data

        merged_bundle = {
            "resourceType": "Bundle",
            "id": str(uuid.uuid4()),
            "meta": {
                "lastUpdated": datetime.now(timezone.utc).isoformat(),
                "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Bundle"]
            },
            "type": "transaction", 
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "entry": []
        }

        org_resource = create_organization_resource(org_data)
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": org_resource,
            "request": {
                "method": "PUT",
                "url": f"Organization/{org_resource['id']}"
            }
        })

        practitioner_resource = create_practitioner_resource(practitioner_data)
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": practitioner_resource,
            "request": {
                "method": "PUT",
                "url": f"Practitioner/{practitioner_resource['id']}"
            }
        })

        role_resource = create_practitioner_role_resource(practitioner_data, org_data)
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": role_resource,
            "request": {
                "method": "PUT",
                "url": f"PractitionerRole/{role_resource['id']}"
            }
        })

        for sample_id, sample_clinical_data in matched_samples.items():
            debug_print(f"Adding patient for sample: {sample_id}")

            patient_resource = create_patient_resource(sample_id, sample_clinical_data, org_data)
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": patient_resource,
                "request": {
                    "method": "PUT",
                    "url": f"Patient/{patient_resource['id']}"
                }
            })

            specimen_resource = create_specimen_resource(sample_id, sample_clinical_data, practitioner_data, org_data)
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": specimen_resource,
                "request": {
                    "method": "PUT",
                    "url": f"Specimen/{specimen_resource['id']}"
                }
            })

            service_request_resource = create_service_request_resource(sample_id, sample_clinical_data, practitioner_data, org_data)
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": service_request_resource,
                "request": {
                    "method": "PUT",
                    "url": f"ServiceRequest/{service_request_resource['id']}"
                }
            })

        observations_by_sample = {}
        for obs in all_observations:
            subject_ref = obs.get('subject', {}).get('reference', '')
            if subject_ref.startswith('Patient/'):
                sample_id = subject_ref.replace('Patient/', '').replace('-patient', '')
                if sample_id not in observations_by_sample:
                    observations_by_sample[sample_id] = []
                observations_by_sample[sample_id].append(obs)

        for sample_id, sample_observations in observations_by_sample.items():
            sample_clinical_data = matched_samples.get(sample_id)
            
            diagnostic_report = create_diagnostic_report(
                sample_id,
                sample_observations,
                sample_clinical_data,
                org_data,
                practitioner_data
            )
            
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": diagnostic_report,
                "request": {
                    "method": "PUT",
                    "url": f"DiagnosticReport/{diagnostic_report['id']}"
                }
            })

        for entry in fhir_bundle.get('entry', []):
            resource = entry.get('resource', {})
            resource_type = resource.get('resourceType')
            resource_id = resource.get('id')
            
            entry_with_request = {
                "fullUrl": entry.get('fullUrl', f"urn:uuid:{str(uuid.uuid4())}"),
                "resource": resource,
                "request": {
                    "method": "PUT",
                    "url": f"{resource_type}/{resource_id}" if resource_id else f"{resource_type}"
                }
            }
            merged_bundle['entry'].append(entry_with_request)

        with open(args.output, 'w') as f:
            json.dump(merged_bundle, f, indent=2)

    except Exception as e:
        debug_print(f"Error occurred: {str(e)}")
        import traceback
        debug_print(f"Traceback: {traceback.format_exc()}")
        sys.exit(1)

if __name__ == "__main__":
    main()
