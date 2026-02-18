#!/usr/bin/env python3

import json
import argparse
import uuid
import os
import sys
import re
from collections import defaultdict
from datetime import datetime, timezone

try:
    import pandas as pd
    import numpy as np
except Exception as e:
    raise SystemExit("pandas required (pip install pandas openpyxl)")

def get_ncbi_gene_id(gene_name):
    mapping = {
        'katG': '885638',
        'rpoB': '888164',
        'inhA': '886364',
        'gyrA': '887103',
        'gyrB': '887081',
        'embB': '886404',
        'pncA': '886105',
        'rpsL': '887363',
        'rrs': '887376',
        'ethA': '885946',
        'ahpC': '888252',
        'fabG1': '886362',
        'eis': '886915',
        'tlyA': '886311',
        'alr': '886443',
        'ddl': '886444',
        'rpoC': '888165',
        'embA': '886405',
        'embC': '886403',
        'ubiA': '887985',
        'ndh': '886523',
        'Rv0678': '888235',
        'gid': '886243'
    }
    return mapping.get(str(gene_name).strip())

def get_drug_panel_config():
    return [
        ('rifampicin', '89489-9', 'rifAMPin [Susceptibility] by Genotype method'),
        ('isoniazid', '89488-1', 'Isoniazid [Susceptibility] by Genotype method'),
        ('ethambutol', '89491-5', 'Ethambutol [Susceptibility] by Genotype method'),
        ('pyrazinamide', '92242-7', 'Pyrazinamide [Susceptibility] by Genotype method'),
        ('moxifloxacin', '96112-8', 'Moxifloxacin [Susceptibility] by Genotype method'),
        ('levofloxacin', '20629-2', 'levoFLOXacin [Susceptibility]'),
        ('bedaquiline', '96107-8', 'Bedaquiline [Susceptibility] by Genotype method'),
        ('delamanid', '96109-4', 'Delamanid [Susceptibility] by Genotype method'),
        ('pretomanid', '93850-6', 'Pretomanid [Susceptibility]'),
        ('streptomycin', '96114-4', 'Streptomycin [Susceptibility] by Genotype method'),
        ('amikacin', '89484-0', 'Amikacin [Susceptibility] by Genotype method'),
        ('kanamycin', '89482-4', 'Kanamycin [Susceptibility] by Genotype method'),
        ('capreomycin', '89483-2', 'Capreomycin [Susceptibility] by Genotype method'),
        ('clofazimine', '96108-6', 'Clofazimine [Susceptibility] by Genotype method'),
        ('ethionamide', '96110-2', 'Ethionamide [Susceptibility] by Genotype method'),
        ('linezolid', '96111-0', 'Linezolid [Susceptibility] by Genotype method'),
        ('cycloserine', '103959-3', 'cycloSERINE [Susceptibility] by Genotype method')
    ]

def normalize_drug_name(drug_name):
    d = str(drug_name).strip().upper()
    if 'RIF' in d: return 'rifampicin'
    if 'INH' in d: return 'isoniazid'
    if 'PZA' in d: return 'pyrazinamide'
    if 'EMB' in d: return 'ethambutol'
    if 'STR' in d: return 'streptomycin'
    if 'LFX' in d: return 'levofloxacin'
    if 'MFX' in d: return 'moxifloxacin'
    if 'KAN' in d: return 'kanamycin'
    if 'AMK' in d: return 'amikacin'
    if 'CAP' in d: return 'capreomycin'
    if 'ETH' in d: return 'ethionamide'
    if 'LIN' in d: return 'linezolid'
    if 'BDQ' in d: return 'bedaquiline'
    if 'CFZ' in d: return 'clofazimine'
    if 'FQ' in d: return 'fluoroquinolone'
    return d.lower()

def clean_variant_text(text):
    if not text: return ""
    cleaned = re.sub(r'\s*\n?\(.*?\)', '', str(text)).strip()
    return cleaned

def get_deeplex_classification_coding(confidence_text):
    if not confidence_text:
        return None
        
    text = str(confidence_text).strip()
    text_lower = text.lower()
    
    if "uncertain" in text_lower:
        return {
            'system': 'http://loinc.org',
            'code': 'LA26333-7',
            'display': 'Uncertain significance'
        }
    
    if "not associated with resistance" in text_lower:
        if "interim" in text_lower:
             return {
                'system': 'http://terminology.kemkes.go.id/sp',
                'code': 'SP000480',
                'display': 'Not assoc w R - Interim'
            }
        return {
            'system': 'http://terminology.kemkes.go.id/sp',
            'code': 'SP000481',
            'display': 'Not assoc w R'
        }
        
    if "associated with resistance" in text_lower:
        if "interim" in text_lower:
             return {
                'system': 'http://terminology.kemkes.go.id/sp',
                'code': 'SP000479',
                'display': 'Assoc w R - Interim'
            }
        return {
            'system': 'http://terminology.kemkes.go.id/sp',
            'code': 'SP000478',
            'display': 'Assoc w R'
        }
    
    return None

def create_observation(sample_id, variant_data, obs_index):
    gene = variant_data.get('gene')
    codon_change = clean_variant_text(variant_data.get('codon_change'))
    amino_acid = clean_variant_text(variant_data.get('amino_acid'))
    confidence = variant_data.get('confidence')
    variant_type = variant_data.get('type', 'resistance') 
    
    components = []

    if gene:
        gene_id = get_ncbi_gene_id(gene)
        gene_code = gene_id if gene_id else gene
        
        components.append({
            "code": {"coding": [{"system": "http://loinc.org", "code": "48018-6", "display": "Gene studied [ID]"}]},
            "valueCodeableConcept": {
                "coding": [{
                    "system": "https://www.ncbi.nlm.nih.gov/gene",
                    "code": gene_code,
                    "display": str(gene)
                }],
                "text": str(gene)
            }
        })

    if codon_change:
        display_val = codon_change
        components.append({
            "code": {"coding": [{"system": "http://loinc.org", "code": "81290-9", "display": "Genomic DNA change (gHGVS)"}]},
            "valueCodeableConcept": {
                "coding": [{"system": "https://varnomen.hgvs.org", "code": display_val, "display": display_val}],
                "text": display_val
            }
        })

    if amino_acid:
        display_aa = amino_acid
        if display_aa.strip().lower().startswith('p.'):
            display_aa_clean = display_aa.strip()[2:]
        else:
            display_aa_clean = display_aa

        if display_aa.lower() not in ['nan', 'uncharact.', 'frameshift']:
            components.append({
                "code": {"coding": [{"system": "http://loinc.org", "code": "48005-3", "display": "Amino acid change (pHGVS)"}]},
                "valueCodeableConcept": {
                    "coding": [{"system": "https://varnomen.hgvs.org", 
                                "code": f"NC_000962.3:p.({display_aa_clean})", 
                                "display": display_aa_clean}],
                    "text": display_aa
                }
            })

    classification_coding = get_deeplex_classification_coding(confidence)
    if classification_coding:
        components.append({
            "code": {"coding": [{"system": "http://loinc.org", "code": "53037-8", "display": "Genetic variation clinical significance [Imp]"}]},
            "valueCodeableConcept": {
                "coding": [{"system": classification_coding['system'], "code": classification_coding['code'], "display": classification_coding['display']}],
                "text": str(confidence)
            }
        })
    
    obs_id = f"{sample_id}-obs-{obs_index}"
    effective_datetime = datetime.now(timezone.utc).isoformat()
    
    return {
        "resourceType": "Observation",
        "id": obs_id,
        "meta": {
            "profile": ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"],
            "tag": [{"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}]
        },
        "status": "final",
        "category": [
            {
                "coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category", "code": "laboratory", "display": "Laboratory"}]
            },
            {
                "coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "GE", "display": "Genetics"}]
            }
        ],
        "code": {"coding": [{"system": "http://loinc.org", "code": "69548-6", "display": "Genetic variant assessment"}]},
        "subject": {"reference": f"Patient/{sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{sample_id}-specimen"},
        "effectiveDateTime": effective_datetime,
        "performer": [{"reference": "Organization/100007732"}],
        "valueCodeableConcept": {"coding": [{"system": "http://loinc.org", "code": "LA9633-4", "display": "Present"}], "text": "Present"},
        "component": components
    }

def create_lineage_observation(sample_id, lineage_text, obs_index):
    obs_id = f"{sample_id}-lineage" 
    effective_datetime = datetime.now(timezone.utc).isoformat()
    
    return {
        "resourceType": "Observation",
        "id": obs_id,
        "status": "final",
        "meta": {
             "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"], 
             "tag": [{"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}]
        },
        "category": [
            {
                "coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category", "code": "laboratory", "display": "Laboratory"}]
            },
            {
                "coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "GE", "display": "Genetics"}]
            }
        ],
        "code": {"coding": [{"system": "http://loinc.org", "code": "614-8", "display": "Mycobacterial strain [Type] in Isolate by Mycobacterial subtyping"}]},
        "subject": {"reference": f"Patient/{sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{sample_id}-specimen"},
        "effectiveDateTime": effective_datetime,
        "performer": [{"reference": "Organization/100007732"}],
        "valueCodeableConcept": {
            "coding": [
                {
                    "system": "http://tb-lineage.org",
                    "code": str(lineage_text),
                    "display": f"TB Lineage {str(lineage_text)}"
                }
            ],
            "text": f"Lineage {str(lineage_text)}"
        }
    }

def create_drug_panel_observation(sample_id, resistant_drugs_detected):
    panel_components = []
    panel_config = get_drug_panel_config()
    
    for drug_key, code, display in panel_config:
        is_resistant = drug_key in resistant_drugs_detected
        
        value_code = "LA6676-6" if is_resistant else "LA24225-7"
        value_display = "Resistant" if is_resistant else "Susceptible"
        
        panel_components.append({
            "code": {
                "coding": [{
                    "system": "http://loinc.org",
                    "code": code,
                    "display": display
                }]
            },
            "valueCodeableConcept": {
                "coding": [{
                    "system": "http://loinc.org",
                    "code": value_code,
                    "display": value_display
                }]
            }
        })

    return {
        "resourceType": "Observation",
        "id": f"{sample_id}-susceptibility-panel",
        "meta": {
            "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"],
        },
        "text": {
            "status": "generated",
            "div": f"<div xmlns=\"http://www.w3.org/1999/xhtml\">Mycobacterial susceptibility panel for {sample_id}</div>"
        },
        "status": "final",
        "category": [
            {
                "coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category", "code": "laboratory", "display": "Laboratory"}]
            },
            {
                "coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074", "code": "GE", "display": "Genetics"}]
            }
        ],
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": "89486-5",
                "display": "Mycobacterial susceptibility panel Qualitative by Genotype method"
            }]
        },
        "subject": {"reference": f"Patient/{sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{sample_id}-specimen"},
        "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
        "performer": [{"reference": "Organization/100007732"}],
        "component": panel_components
    }

def process_deeplex_batch(input_path, output_dir):
    all_sample_data = defaultdict(lambda: {"display_id": None, "variants": [], "lineages": set(), "resistant_drugs": set()})
    
    sheets_to_process = [
        {'name': 'Drug resistance variants', 'type': 'resistance'},
        {'name': 'Uncharacterised variants', 'type': 'uncharacterised'},
        {'name': 'Phylo_Syn variants', 'type': 'phylo'}
    ]

    try:
        xl = pd.ExcelFile(input_path)
        
        for sheet_info in sheets_to_process:
            sheet_name = sheet_info['name']
            var_type = sheet_info['type']
            
            if sheet_name not in xl.sheet_names:
                print(f"Sheet '{sheet_name}' not found. Skipping.")
                continue
                
            print(f"Processing sheet: {sheet_name}")
            df = pd.read_excel(input_path, sheet_name=sheet_name, header=None)
            
            header_drug = df.iloc[0].fillna(method='ffill')
            header_gene = df.iloc[1].fillna(method='ffill')
            header_col = df.iloc[2]
            
            lineage_col_idx = -1
            if var_type == 'resistance':
                for i, col_name in enumerate(header_drug):
                    if "SNP-based lineage" in str(col_name) or "SNP-based lineage" in str(header_col[i]):
                        lineage_col_idx = i
                        break

            for idx, row in df.iterrows():
                if idx < 3: continue
                
                raw_id = row[0]
                if pd.isna(raw_id): continue
                
                sample_id = str(raw_id).strip()
                if isinstance(raw_id, float) and raw_id.is_integer():
                    sample_id = str(int(raw_id))
                elif sample_id.endswith('.0'):
                    sample_id = sample_id[:-2]
                
                if not sample_id or sample_id.lower() == 'nan': continue
                
                norm_id = sample_id.lower()
                norm_id = re.sub(r'^sample\s*', '', norm_id)
                core_id = re.sub(r'[^a-z0-9]', '', norm_id)
                
                if not core_id: continue 

                current_display = all_sample_data[core_id]["display_id"]
                if current_display is None:
                    all_sample_data[core_id]["display_id"] = sample_id
                else:
                    if sample_id.lower().startswith('sample') and not current_display.lower().startswith('sample'):
                        all_sample_data[core_id]["display_id"] = sample_id

                for i in range(len(header_col)):
                    col_name = str(header_col[i]).strip().lower()
                    
                    if "codon change" in col_name:
                        val_codon_raw = str(row[i]).strip()
                        
                        if not val_codon_raw or val_codon_raw.lower() in ['nan', 'wt', 'wt*', 'low coverage', 'nd', 'susceptible']: 
                            continue
                        
                        val_aa_raw = ""
                        val_conf_raw = ""
                        val_res_raw = ""
                        
                        for offset in range(1, 5):
                            if i + offset >= len(header_col): break
                            next_header = str(header_col[i+offset]).strip().lower()
                            
                            if "amino acid" in next_header:
                                val_aa_raw = str(row[i+offset]).strip()
                            elif "confidence" in next_header:
                                val_conf_raw = str(row[i+offset]).strip()
                            elif "resistance" in next_header:
                                val_res_raw = str(row[i+offset]).strip()
                        
                        codons = val_codon_raw.split('|')
                        aas = val_aa_raw.split('|')
                        confs = val_conf_raw.split('|')
                        ress = val_res_raw.split('|')
                        
                        count = len(codons)
                        
                        def get_val(lst, idx):
                            if not lst: return ""
                            if idx < len(lst): return lst[idx].strip()
                            if len(lst) > 0: return lst[-1].strip()
                            return ""

                        for k in range(count):
                            c_change = codons[k].strip()
                            if not c_change or c_change.lower() in ['nan', 'wt', 'wt*', 'low coverage', 'nd', 'susceptible']:
                                continue
                            
                            conf_val = get_val(confs, k)
                            drug_name = str(header_drug[i]).strip()
                            normalized_drug = normalize_drug_name(drug_name)

                            if var_type == 'resistance' and "associated with resistance" in str(conf_val).lower() and "not associated" not in str(conf_val).lower():
                                all_sample_data[core_id]["resistant_drugs"].add(normalized_drug)
                                if 'fluoroquinolone' in normalized_drug:
                                    all_sample_data[core_id]["resistant_drugs"].add('moxifloxacin')
                                    all_sample_data[core_id]["resistant_drugs"].add('levofloxacin')
                                    all_sample_data[core_id]["resistant_drugs"].add('ciprofloxacin')
                                    all_sample_data[core_id]["resistant_drugs"].add('ofloxacin')
                                
                            variant_data = {
                                'drug_raw': drug_name,
                                'gene': str(header_gene[i]).strip(),
                                'codon_change': c_change,
                                'amino_acid': get_val(aas, k),
                                'confidence': conf_val,
                                'resistance': get_val(ress, k),
                                'type': var_type
                            }
                            all_sample_data[core_id]["variants"].append(variant_data)

                if lineage_col_idx != -1:
                    lineage_val = str(row[lineage_col_idx]).strip()
                    if lineage_val and lineage_val.lower() != 'nan':
                        all_sample_data[core_id]["lineages"].add(lineage_val)

        processed_count = 0
        print("\n--- Processing Summary ---")
        for core_id, data in all_sample_data.items():
            final_display_id = data["display_id"]
            safe_id = re.sub(r'[^a-z0-9_-]', '', final_display_id.lower())
            
            observations = []
            obs_counter = 1
            
            for var_data in data["variants"]:
                observations.append(create_observation(safe_id, var_data, obs_counter))
                obs_counter += 1
            
            observations.append(create_drug_panel_observation(safe_id, data["resistant_drugs"]))
            
            for lin_text in data["lineages"]:
                observations.append(create_lineage_observation(safe_id, lin_text, obs_counter))
                obs_counter += 1
            
            bundle = {
                "resourceType": "Bundle",
                "type": "collection",
                "entry": [{"fullUrl": f"urn:uuid:{obs['id']}", "resource": obs} for obs in observations]
            }
            
            out_file = os.path.join(output_dir, f"{safe_id}.fhir.json")
            with open(out_file, 'w', encoding='utf-8') as f:
                json.dump(bundle, f, indent=2, ensure_ascii=False)
            
            processed_count += 1
            print(f"Sample: {final_display_id} (Merged ID: {core_id}) -> {len(observations)} observations")

        print(f"\nSuccessfully processed {processed_count} unique samples.")

    except Exception as e:
        print(f"Error processing batch: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Deeplex XLSX file')
    parser.add_argument('--output_dir', required=True, help='Output directory for JSON files')
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    process_deeplex_batch(args.input, args.output_dir)

if __name__ == "__main__":
    main()