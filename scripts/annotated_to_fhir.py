#!/usr/bin/env python3

import json
import csv
import argparse
from datetime import datetime, timezone
import os
import glob
from collections import defaultdict
import sys
import uuid
import re

def load_lineage_data(lineage_dir):
    lineage_data = {}
    
    if not os.path.exists(lineage_dir):
        return lineage_data
    
    lineage_files = glob.glob(os.path.join(lineage_dir, "*.lineage.json"))
    
    for lineage_file in lineage_files:
        try:
            with open(lineage_file, 'r') as f:
                lineage_info = json.load(f)
                
                sample_id = lineage_info.get('sample_id', 
                    os.path.basename(lineage_file).replace('.lineage.json', ''))
                
                possible_ids = [
                    sample_id,
                    sample_id.replace('_ont.fastq', ''),
                    sample_id.replace('_illumina', ''),
                    sample_id.replace('.fastq', ''),
                    sample_id.split('_')[0], 
                    os.path.basename(lineage_file).replace('.lineage.json', '')
                ]
                
                for pid in set(possible_ids):
                    if pid:
                        lineage_data[pid] = lineage_info
                
        except Exception as e:
            print(f"Error loading lineage file {lineage_file}: {e}")
            continue
    
    return lineage_data

def fix_malformed_vcf(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line_num, line in enumerate(infile, 1):
            if line.startswith('#'):
                outfile.write(line)
            else:
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    if len(fields) >= 9 and ('significance' in fields[8] or 'VARIANT_ID' in fields[8]):
                        outfile.write('\t'.join(fields[:8]) + '\n')
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)

def simple_vcf_parser(vcf_file):
    class SimpleVCFRecord:
        def __init__(self, chrom, pos, ref, alt, info_dict):
            self.CHROM = chrom
            self.POS = int(pos)
            self.REF = ref
            self.ALT = [alt]
            self.INFO = info_dict
    
    variants = []
    
    with open(vcf_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) >= 8:
                chrom = fields[0]
                pos = fields[1]
                ref = fields[3]
                alt = fields[4]
                info_str = fields[7]
                
                # Parse INFO field
                info_dict = {}
                if info_str != '.':
                    for info_item in info_str.split(';'):
                        if '=' in info_item:
                            key, value = info_item.split('=', 1)
                            info_dict[key] = value
                        else:
                            info_dict[info_item] = True
                
                variants.append(SimpleVCFRecord(chrom, pos, ref, alt, info_dict))
    
    return variants

def get_variant_value(record, field, default='unknown'):
    try:
        if field in record.INFO:
            value = record.INFO[field]
            return value if value is not None else default
        return default
    except Exception as e:
        return default

def sanitize_id(id_string):
    sanitized = re.sub(r'[^A-Za-z0-9\-]', '-', str(id_string))
    sanitized = re.sub(r'-+', '-', sanitized)
    sanitized = sanitized.strip('-')
    return sanitized if sanitized else "unknown"

def extract_sample_id_from_filename(filename):
    basename = os.path.basename(filename)
    
    basename = basename.replace('.annotated_variants.vcf.gz', '')
    basename = basename.replace('.annotated_variants.vcf', '')
    basename = basename.replace('.vcf.gz', '')
    basename = basename.replace('.vcf', '')
    
    clean_basename = basename.replace('_ont', '').replace('_illumina', '').replace('.fastq', '')
    if '_' in clean_basename:
        clean_basename = clean_basename.split('_')[0]
    
    variations = [
        clean_basename, 
        basename,
        basename.replace('_ont', ''),
        basename.replace('_illumina', ''),
        basename.replace('.fastq', ''),
        basename.split('_')[0] if '_' in basename else basename 
    ]
    
    return variations

def get_drug_snomed_mapping(drug_name):
    drug_mapping = {
        'rifampicin': {'code': '29175007', 'display': 'Product containing rifampicin'},
        'isoniazid': {'code': '81335000', 'display': 'Product containing isoniazid'},
        'pyrazinamide': {'code': '13592004', 'display': 'Product containing pyrazinamide'},
        'ethambutol': {'code': '24450004', 'display': 'Product containing ethambutol'},
        'streptomycin': {'code': '40877002', 'display': 'Product containing streptomycin'},
        'fluoroquinolone': {'code': '1010205001', 'display': 'Medicinal product containing fluoroquinolone and acting as antibacterial agent'},
        'levofloxacin': {'code': '96087006', 'display': 'Product containing levofloxacin'},
        'moxifloxacin': {'code': '371296007', 'display': 'Product containing moxifloxacin'},
        'ofloxacin': {'code': '96086002', 'display': 'Product containing ofloxacin'},
        'ciprofloxacin': {'code': '7577004', 'display': 'Product containing ciprofloxacin'},
        'gatifloxacin': {'code': '371238005', 'display': 'Product containing gatifloxacin'},
        'amikacin': {'code': '48836000', 'display': 'Product containing amikacin'},
        'kanamycin': {'code': '71451001', 'display': 'Product containing kanamycin'},
        'capreomycin': {'code': '14170004', 'display': 'Product containing capreomycin'},
        'ethionamide': {'code': '414148003', 'display': 'Product containing ethionamide'},
        'linezolid': {'code': '125695009', 'display': 'Product containing linezolid'},
        'bedaquiline': {'code': '714087008', 'display': 'Product containing bedaquiline'},
        'clofazimine': {'code': '72924009', 'display': 'Product containing clofazimine'}
    }
    
    normalized_drug = drug_name.lower().strip()
    return drug_mapping.get(normalized_drug)

def get_who_classification_coding(who_classification):
    classification_mapping = {
        'Uncertain significance': {
            'system': 'http://loinc.org',
            'code': 'LA26333-7',
            'display': 'Uncertain significance'
        },
        'Assoc w R': {
            'system': 'http://terminology.kemkes.go.id/sp',
            'code': 'SP000478',
            'display': 'Assoc w R'
        },
        'Assoc w R - Interim': {
            'system': 'http://terminology.kemkes.go.id/sp',
            'code': 'SP000479',
            'display': 'Assoc w R - Interim'
        },
        'Not assoc w R - Interim': {
            'system': 'http://terminology.kemkes.go.id/sp',
            'code': 'SP000480',
            'display': 'Not assoc w R - Interim'
        },
        'Not assoc w R': {
            'system': 'http://terminology.kemkes.go.id/sp',
            'code': 'SP000481',
            'display': 'Not assoc w R'
        }
    }
    
    return classification_mapping.get(who_classification.strip())

def get_effect_so_mapping(effect):
    effect_mapping = {
        'missense_variant': {'code': 'SO:0001583', 'display': 'missense_variant'},
        'synonymous_variant': {'code': 'SO:0001819', 'display': 'synonymous_variant'},
        'stop_gained': {'code': 'SO:0001587', 'display': 'stop_gained'},
        'stop_lost': {'code': 'SO:0001578', 'display': 'stop_lost'},
        'frameshift_variant': {'code': 'SO:0001589', 'display': 'frameshift_variant'},
        'inframe_insertion': {'code': 'SO:0001821', 'display': 'inframe_insertion'},
        'inframe_deletion': {'code': 'SO:0001822', 'display': 'inframe_deletion'},
        'splice_site_variant': {'code': 'SO:0001629', 'display': 'splice_site_variant'},
        'upstream_gene_variant': {'code': 'SO:0001631', 'display': 'upstream_gene_variant'},
        'downstream_gene_variant': {'code': 'SO:0001632', 'display': 'downstream_gene_variant'},
        'intergenic_variant': {'code': 'SO:0001628', 'display': 'intergenic_variant'},
        'intron_variant': {'code': 'SO:0001627', 'display': 'intron_variant'},
        '5_prime_UTR_variant': {'code': 'SO:0001623', 'display': '5_prime_UTR_variant'},
        '3_prime_UTR_variant': {'code': 'SO:0001624', 'display': '3_prime_UTR_variant'},
        'start_lost': {'code': 'SO:0002012', 'display': 'start_lost'},
        'stop_retained_variant': {'code': 'SO:0001567', 'display': 'stop_retained_variant'},
        'protein_altering_variant': {'code': 'SO:0001818', 'display': 'protein_altering_variant'},
        'coding_sequence_variant': {'code': 'SO:0001580', 'display': 'coding_sequence_variant'},
        'non_coding_transcript_variant': {'code': 'SO:0001619', 'display': 'non_coding_transcript_variant'},
        'regulatory_region_variant': {'code': 'SO:0001566', 'display': 'regulatory_region_variant'}
    }
    
    normalized_effect = effect.lower().strip()
    return effect_mapping.get(normalized_effect)

def get_ncbi_gene_id(gene_name):
    mapping = {
        'katG': '885638',
        'rpoB': '888164',
        'inhA': '886523',
        'gyrA': '887105',
        'gyrB': '887081',
        'embB': '886126',
        'pncA': '888260',
        'rpsL': '888259',
        'rrs': '2700429',
        'ethA': '886175',
        'ahpC': '885717',
        'fabG1': '886551',
        'eis': '885903',
        'tlyA': '885396',
        'alr': '887634',
        'ddlA': '888415',
        'rpoC': '888177',
        'embA': '886123',
        'embC': '886112',
        'ubiA': '886129',
        'ndh': '885746',
        'ndhA': '886430',
        'Rv0678': '888235',
        'gid': '886243'
    }
    return mapping.get(gene_name)

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

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='Path to VCF file')
parser.add_argument('--output', required=True, help='Path to output FHIR JSON')
parser.add_argument('--lineage_dir', help='Directory containing lineage JSON files')
args = parser.parse_args()

lineage_data = {}
if args.lineage_dir:
    lineage_data = load_lineage_data(args.lineage_dir)

file_sample_id_variations = extract_sample_id_from_filename(args.input)

matched_sample_id = None
sample_lineage_info = None
for variation in file_sample_id_variations:
    if variation in lineage_data:
        matched_sample_id = variation
        sample_lineage_info = lineage_data[variation]
        break

if sample_lineage_info:
    print(f"Found lineage info for sample: {matched_sample_id}")
else:
    print("No lineage info found for this sample")

try:
    if not os.path.exists(args.input):
        sys.exit(1)
    
    if os.path.getsize(args.input) == 0:
        sys.exit(1)

    fixed_vcf = args.input + '.fixed'
    fix_malformed_vcf(args.input, fixed_vcf)

    filename = os.path.basename(args.input).lower()
    
    variants = simple_vcf_parser(fixed_vcf)
    
    bundles = defaultdict(list)
    variant_count = 0
    successful_annotations = 0
    
    resistant_drugs_detected = set()

    for idx, record in enumerate(variants):
        try:
            file_sample_id = extract_sample_id_from_filename(args.input)[0]  
            sample_id = f"{file_sample_id}-patient"  
            pos = str(record.POS)
            ref = record.REF
            alt = str(record.ALT[0])
            
            gene = get_variant_value(record, 'GENE')
            drug = get_variant_value(record, 'DRUG')
            effect = get_variant_value(record, 'EFFECT')
            
            who_classification = get_variant_value(record, 'WHO_CLASSIFICATION')
            
            variant_id = get_variant_value(record, 'VARIANT_ID')
            genome_position = get_variant_value(record, 'GENOME_POSITION')
            depth_val = get_variant_value(record, 'DP', None)

            if drug != "unknown" and who_classification != "unknown":
                if "Assoc w R" in who_classification:
                    d_lower = drug.lower().strip()
                    resistant_drugs_detected.add(d_lower)
                    if 'fluoroquinolone' in d_lower:
                        resistant_drugs_detected.add('moxifloxacin')
                        resistant_drugs_detected.add('levofloxacin')
                        resistant_drugs_detected.add('ciprofloxacin')
                        resistant_drugs_detected.add('ofloxacin')

            if gene != "unknown" or drug != "unknown" or who_classification != "unknown":
                successful_annotations += 1

            components = [
                {
                    "code": { 
                      "coding": [
                        {
                          "system": "http://loinc.org",
                          "code": "81290-9",
                          "display": "Genomic DNA change (gHGVS)"
                        }
                      ]
                    },
                    "valueCodeableConcept": {
                        "coding": [{
                            "system": "https://varnomen.hgvs.org",
                            "code": f"NC_000962.3:g.{pos}{ref}>{alt}",
                            "display": f"NC_000962.3:g.{pos}{ref}>{alt}"
                        }]
                    }
                },
                {
                    "code": {"coding": [{"system": "http://www.ncbi.nlm.nih.gov/refseq", "code": "NC_000962.3"}]}
                }
            ]

            if depth_val is not None and depth_val != 'unknown':
                try:
                    components.append({
                        "code": {"coding": [{"system": "http://loinc.org", "code": "82121-5", "display": "Allelic read depth"}]},
                        "valueQuantity": {
                            "value": int(depth_val),
                            "unit": "reads per base pair",
                            "system": "http://unitsofmeasure.org",
                            "code": "[1]"
                        }
                    })
                except ValueError:
                    pass

            if gene != "unknown":
                gene_id = get_ncbi_gene_id(gene)
                
                gene_code = gene_id if gene_id else gene
                
                components.append({
                    "code": {"coding": [{"system": "http://loinc.org", "code": "48018-6", "display": "Gene studied [ID]"}]},
                    "valueCodeableConcept": {
                        "coding": [{
                            "system": "https://www.ncbi.nlm.nih.gov/gene",
                            "code": gene_code,
                            "display": gene
                        }],
                        "text": gene
                    }
                })

            if effect != "unknown":
                effect_mapping = get_effect_so_mapping(effect)
                if effect_mapping:
                    components.append({
                        "code": {"coding": [{"system": "http://loinc.org", "code": "48019-4", "display": "DNA change type"}]},
                        "valueCodeableConcept": {
                            "coding": [{
                                "system": "http://www.sequenceontology.org",
                                "code": effect_mapping['code'],
                                "display": effect_mapping['display']
                            }],
                            "text": effect
                        }
                    })
                else:
                    components.append({
                        "code": {"coding": [{"system": "http://loinc.org", "code": "48019-4", "display": "DNA change type"}]},
                        "valueCodeableConcept": {
                            "text": effect
                        }
                    })

            if who_classification != "unknown":
                who_coding = get_who_classification_coding(who_classification)
                if who_coding:
                    components.append({
                        "code": {"coding": [{"system": "http://loinc.org", "code": "53037-8", "display": "Genetic variation clinical significance [Imp]"}]},
                        "valueCodeableConcept": {
                            "coding": [{
                                "system": who_coding['system'],
                                "code": who_coding['code'],
                                "display": who_coding['display']
                            }],
                            "text": who_classification
                        }
                    })
                else:
                    components.append({
                        "code": {"coding": [{"system": "http://loinc.org", "code": "53037-8", "display": "Genetic variation clinical significance [Imp]"}]},
                        "valueCodeableConcept": {
                            "text": who_classification
                        }
                    })

            if variant_id != "unknown":
                clean_variant_id = variant_id.split('_')[-1] if '_' in variant_id else variant_id
                
                if clean_variant_id.startswith('p.'):
                    clean_variant_id = clean_variant_id[2:]

                components.append({
                    "code": {"coding": [{"system": "http://loinc.org", "code": "48005-3", "display": "Amino acid change (pHGVS)"}]},
                    "valueCodeableConcept": {
                        "coding": [{
                            "system": "https://varnomen.hgvs.org",
                            "code": f"NC_000962.3:p.({clean_variant_id})",
                            "display": clean_variant_id
                        }]
                    }
                })

            try:
                position_value = int(pos)
                components.append({
                    "code": {"coding": [{"system": "http://loinc.org", "code": "81254-5", "display": "Variant exact start-end"}]},
                    "valueRange": {
                        "low": {
                            "value": position_value
                        }
                    }
                })
            except ValueError:
                pass

            annotation_text = f"Genomic variant at position {pos}: {ref}>{alt}"
            if gene != "unknown":
                annotation_text += f" in gene {gene}"
            if effect != "unknown":
                annotation_text += f" ({effect})"
            if variant_id != "unknown":
                annotation_text += f" - {variant_id}"
            if drug != "unknown":
                 annotation_text += f" - Associated with {drug}"
            if who_classification != "unknown":
                annotation_text += f" - WHO Classification: {who_classification}"
            
            div_text = f"<div xmlns=\"http://www.w3.org/1999/xhtml\">{annotation_text}</div>"

            observation_id = f"{file_sample_id}-obs-{idx+1}"
            observation = {
                "resourceType": "Observation",
                "id": observation_id,
                "meta": {
                    "profile": [
                        "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/variant"
                    ],
                    "tag": [
                        {
                            "system": "http://terminology.kemkes.go.id/sp",
                            "code": "genomics",
                            "display": "Genomics"
                        }
                    ]
                },
                "text": {
                    "status": "generated",
                    "div": div_text
                },
                "status": "final",
                "category": [
                    {
                        "coding": [{
                            "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                            "code": "laboratory",
                            "display": "Laboratory"
                        }]
                    },
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
                        "code": "69548-6",
                        "display": "Genetic variant assessment"
                    }]
                },
                "valueCodeableConcept": {
                    "coding": [{
                        "system": "http://loinc.org",
                        "code": "LA9633-4",
                        "display": "Present"
                    }],
                    "text": "Present"
                },
                "subject": {"reference": f"Patient/{file_sample_id}-patient"},
                "specimen": {"reference": f"Specimen/{file_sample_id}-specimen"},
                "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
                "performer": [{"reference": "Organization/100007732"}],
                "component": components
            }
            
            bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": observation
            })

            variant_count += 1
            if variant_count % 100 == 0:
                print(f"Processed {variant_count} variants")
                
        except Exception as e:
            print(f"Error processing variant {idx}: {e}")
            continue
    
    primary_sample_id = extract_sample_id_from_filename(args.input)[0]
    
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

    panel_observation = {
        "resourceType": "Observation",
        "id": f"{primary_sample_id}-susceptibility-panel",
        "meta": {
            "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"],
        },
        "text": {
            "status": "generated",
            "div": f"<div xmlns=\"http://www.w3.org/1999/xhtml\">Mycobacterial susceptibility panel for {primary_sample_id}</div>"
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
        "subject": {"reference": f"Patient/{primary_sample_id}-patient"},
        "specimen": {"reference": f"Specimen/{primary_sample_id}-specimen"},
        "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
        "performer": [{"reference": "Organization/100007732"}],
        "component": panel_components
    }

    bundles[primary_sample_id].append({
        "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
        "resource": panel_observation
    })

    if sample_lineage_info:
        lineage = sample_lineage_info.get('lineage', 'unknown')
        family = sample_lineage_info.get('family', 'Unknown')
        
        if lineage != 'unknown':
            lineage_obs_id = f"{primary_sample_id}-lineage"
            lineage_text = f"Mycobacterial Lineage: {lineage}"
            if family != 'Unknown':
                lineage_text += f" ({family})"
            
            div_text = f"<div xmlns=\"http://www.w3.org/1999/xhtml\">{lineage_text}</div>"

            lineage_observation = {
                "resourceType": "Observation",
                "id": lineage_obs_id,
                "meta": {
                    "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"],
                    "tag": [{"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}]
                },
                "text": {
                    "status": "generated",
                    "div": div_text
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
                        "code": "614-8",
                        "display": "Mycobacterial strain [Type] in Isolate by Mycobacterial subtyping"
                    }]
                },
                "valueCodeableConcept": {
                    "coding": [{
                        "system": "http://tb-lineage.org",
                        "code": lineage,
                        "display": f"TB Lineage {lineage}"
                    }],
                    "text": f"Lineage {lineage}"
                },
                "subject": {"reference": f"Patient/{primary_sample_id}-patient"},
                "specimen": {"reference": f"Specimen/{primary_sample_id}-specimen"},
                "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
                "performer": [{"reference": "Organization/100007732"}]
            }
            
            bundles[primary_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": lineage_observation
            })

    fhir_output = {
        "resourceType": "Bundle",
        "type": "collection",
        "entry": []
    }

    for sample_id, entries in bundles.items():
        fhir_output['entry'].extend(entries)

    with open(args.output, 'w') as out:
        json.dump(fhir_output, out, indent=2)

    try:
        os.remove(fixed_vcf)
    except:
        pass

except Exception as e:
    import traceback
    print(f"Error: {e}")
    traceback.print_exc()
    sys.exit(1)