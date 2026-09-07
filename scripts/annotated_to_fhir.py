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
from clinical_metadata_parser import load_organization_metadata


class AnnotationError(Exception):
    pass


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


def load_coverage_data(coverage_dir):
    coverage_data = {}

    if not coverage_dir or not os.path.exists(coverage_dir):
        return coverage_data

    for coverage_file in glob.glob(os.path.join(coverage_dir, "*.coverage.json")):
        try:
            with open(coverage_file, 'r') as f:
                info = json.load(f)

            sample_id = info.get('sample_id',
                os.path.basename(coverage_file).replace('.coverage.json', ''))

            possible_ids = [
                sample_id,
                sample_id.replace('_ont.fastq', ''),
                sample_id.replace('_illumina', ''),
                sample_id.replace('.fastq', ''),
                sample_id.split('_')[0],
                os.path.basename(coverage_file).replace('.coverage.json', '')
            ]

            for pid in set(possible_ids):
                if pid:
                    coverage_data[pid] = info

        except Exception as e:
            print(f"Error loading coverage file {coverage_file}: {e}")
            continue

    return coverage_data


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
        def __init__(self, chrom, pos, ref, alt, qual, info_dict):
            self.CHROM = chrom
            self.POS = int(pos)
            self.REF = ref
            self.ALT = [alt]
            self.QUAL = qual
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
                qual_str = fields[5]
                info_str = fields[7]

                try:
                    qual = float(qual_str) if qual_str not in ('.', '') else None
                except ValueError:
                    qual = None

                info_dict = {}
                if info_str != '.':
                    for info_item in info_str.split(';'):
                        if '=' in info_item:
                            key, value = info_item.split('=', 1)
                            info_dict[key] = value
                        else:
                            info_dict[info_item] = True

                variants.append(SimpleVCFRecord(chrom, pos, ref, alt, qual, info_dict))

    return variants


def get_variant_value(record, field, default='unknown'):
    try:
        if field in record.INFO:
            value = record.INFO[field]
            return value if value is not None else default
        return default
    except Exception as e:
        return default


def get_variant_list(record, field):
    raw = get_variant_value(record, field, default=None)
    if raw is None or isinstance(raw, bool):
        return []
    if isinstance(raw, (list, tuple)):
        values = [str(v).strip() for v in raw]
    else:
        values = [v.strip() for v in str(raw).split(',')]
    if not values or all(not v or v.lower() == 'unknown' for v in values):
        return []
    return values


def resolve_variant_drug_pairs(record, allow_legacy=False):

    drugs = get_variant_list(record, 'DRUG')
    grades = get_variant_list(record, 'WHO_CLASSIFICATION')
    variant_ids = get_variant_list(record, 'VARIANT_ID')

    if not drugs:
        return []

    if len(grades) == len(drugs):
        pairs = []
        for idx, drug in enumerate(drugs):
            variant_id = variant_ids[idx] if idx < len(variant_ids) else (
                variant_ids[0] if len(variant_ids) == 1 else None)
            pairs.append((drug, grades[idx], variant_id))
        return pairs

    location = f"{record.CHROM}:{record.POS} {record.REF}>{record.ALT[0]}"

    if not allow_legacy:
        raise AnnotationError(
            f"{location}: {len(drugs)} drugs but {len(grades)} WHO grades "
            f"({drugs} / {grades}).\n"
            f"  The annotation table has one row per (variant, drug) pair, so "
            f"bcftools annotate collapsed the grade and the per-drug association is "
            f"no longer recoverable.\n"
            f"  Fix: rebuild the table with scripts/rebuild_annotation_table.sh, then "
            f"point params.annotation_table / params.annotation_header at the result.\n"
            f"  To proceed anyway with the most severe grade applied to every listed "
            f"drug (this over-calls resistance), rerun with --legacy-annotation-table."
        )

    grade = most_severe_grade(grades) if grades else 'unknown'
    variant_id = variant_ids[0] if variant_ids else None
    return [(drug, grade, variant_id) for drug in drugs]


WHO_GRADE_SEVERITY = {
    'assoc w r': 5,
    'assoc w r - interim': 4,
    'uncertain significance': 3,
    'not assoc w r - interim': 2,
    'not assoc w r': 1,
}


def most_severe_grade(grades):
    if not grades:
        return 'unknown'
    return max(grades, key=lambda g: WHO_GRADE_SEVERITY.get(g.strip().lower(), 0))


def is_resistance_grade(grade):
    return grade.strip().lower().startswith('assoc w r')


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
        'clofazimine': {'code': '72924009', 'display': 'Product containing clofazimine'},
        'delamanid': {'code': '714098009', 'display': 'Product containing delamanid'},
        'pretomanid': {'code': '789321008', 'display': 'Product containing pretomanid'},
        'cycloserine': {'code': '51334008', 'display': 'Product containing cycloserine'},
        'terizidone': {'code': '1395908004', 'display': 'Product containing terizidone'},
        'para-aminosalicylic_acid': {'code': '417238004', 'display': 'Product containing aminosalicylic acid'},
    }

    normalized_drug = drug_name.lower().strip().replace(' ', '_')
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
        'regulatory_region_variant': {'code': 'SO:0001566', 'display': 'regulatory_region_variant'},
        'loss_of_function_variant': {'code': 'SO:0002054', 'display': 'loss_of_function_variant'},
        'non_coding_transcript_exon_variant': {'code': 'SO:0001792', 'display': 'non_coding_transcript_exon_variant'},
        'initiator_codon_variant': {'code': 'SO:0001582', 'display': 'initiator_codon_variant'},
        'feature_ablation': {'code': 'SO:0001879', 'display': 'feature_ablation'},
        'feature_truncation': {'code': 'SO:0001906', 'display': 'feature_truncation'},
        'disruptive_inframe_deletion': {'code': 'SO:0001826', 'display': 'disruptive_inframe_deletion'},
        'disruptive_inframe_insertion': {'code': 'SO:0001824', 'display': 'disruptive_inframe_insertion'},
        'splice_acceptor_variant': {'code': 'SO:0001574', 'display': 'splice_acceptor_variant'},
        'splice_donor_variant': {'code': 'SO:0001575', 'display': 'splice_donor_variant'},
        'splice_region_variant': {'code': 'SO:0001630', 'display': 'splice_region_variant'},
        'transcript_ablation': {'code': 'SO:0001893', 'display': 'transcript_ablation'},
        'rrna_variant': {'code': 'SO:0001637', 'display': 'rRNA_gene_variant'},
    }

    normalized_effect = effect.lower().strip()
    return effect_mapping.get(normalized_effect)


_NCBI_GENE_IDS = {
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
    'gid': '886243',
}

_external_gene_map_loaded = False
_unmapped_genes = set()


def load_external_gene_map(path):

    global _external_gene_map_loaded

    if not path or not os.path.exists(path):
        return

    added = 0
    try:
        with open(path) as handle:
            for line in handle:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 2:
                    continue
                gene, gene_id = parts[0].strip(), parts[1].strip()
                if gene and gene_id.isdigit():
                    _NCBI_GENE_IDS[gene] = gene_id
                    added += 1
    except Exception as e:
        print(f"Warning: could not read gene map {path}: {e}")
        return

    _external_gene_map_loaded = True
    print(f"Loaded {added} gene->NCBI mappings from {path}")


def get_ncbi_gene_id(gene_name):
    return _NCBI_GENE_IDS.get(gene_name)


def build_gene_codeable_concept(gene):

    gene_id = get_ncbi_gene_id(gene)

    if gene_id:
        return {
            "coding": [{
                "system": "https://www.ncbi.nlm.nih.gov/gene",
                "code": gene_id,
                "display": gene
            }],
            "text": gene
        }

    _unmapped_genes.add(gene)
    return {"text": gene}


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


def compute_allele_fraction(record):

    dp4 = get_variant_value(record, 'DP4', default=None)
    if dp4 and dp4 != 'unknown':
        try:
            counts = [int(x) for x in str(dp4).split(',')]
            if len(counts) == 4 and sum(counts) > 0:
                return (counts[2] + counts[3]) / sum(counts)
        except (ValueError, ZeroDivisionError):
            pass

    for field in ('AF', 'AF1'):
        value = get_variant_value(record, field, default=None)
        if value and value != 'unknown':
            try:
                fraction = float(str(value).split(',')[0])
                if 0.0 <= fraction <= 1.0:
                    return fraction
            except ValueError:
                continue

    return None


def build_region_studied_observations(sample_id, coverage_info, org_id, panel_drug_keys):

    observations = []

    if not coverage_info or not coverage_info.get('coverage_assessed'):
        return observations

    genes = coverage_info.get('genes', {})
    min_depth = coverage_info.get('min_depth', 10)
    breadth_key = f'breadth_at_{min_depth}x'

    for gene in sorted(genes):
        info = genes[gene]
        gene_drugs = [d for d in info.get('drugs', []) if d.lower() in panel_drug_keys]
        if not gene_drugs:
            continue

        breadth = info.get(breadth_key, 0.0)
        adequate = info.get('adequate', False)

        components = [
            {
                "code": {"coding": [{"system": "http://loinc.org", "code": "48018-6",
                                     "display": "Gene studied [ID]"}]},
                "valueCodeableConcept": build_gene_codeable_concept(gene)
            },
            {
                "code": {"coding": [{"system": "http://loinc.org", "code": "48013-7",
                                     "display": "Genomic reference sequence ID"}]},
                "valueCodeableConcept": {
                    "coding": [{"system": "http://www.ncbi.nlm.nih.gov/refseq",
                                "code": "NC_000962.3", "display": "NC_000962.3"}],
                    "text": "NC_000962.3"
                }
            },
            {
                "code": {"coding": [{"system": "http://loinc.org", "code": "82121-5",
                                     "display": "Allelic read depth"}]},
                "valueQuantity": {
                    "value": info.get('mean_depth', 0),
                    "unit": "reads per base pair",
                    "system": "http://unitsofmeasure.org",
                    "code": "[1]"
                }
            },
            {
                "code": {
                    "coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                "code": "coverage-breadth",
                                "display": f"Fraction of bases at >= {min_depth}x"}],
                    "text": f"Fraction of bases at >= {min_depth}x"
                },
                "valueQuantity": {
                    "value": round(float(breadth) * 100, 2),
                    "unit": "%",
                    "system": "http://unitsofmeasure.org",
                    "code": "%"
                }
            }
        ]

        status_text = "adequately covered" if adequate else "NOT adequately covered"
        div_text = (f"<div xmlns=\"http://www.w3.org/1999/xhtml\">Region studied: {gene} "
                    f"({', '.join(sorted(gene_drugs))}) - {status_text}: "
                    f"{float(breadth) * 100:.1f}% of bases at >= {min_depth}x, "
                    f"mean depth {info.get('mean_depth', 0)}x</div>")

        observations.append({
            "resourceType": "Observation",
            "id": f"{sample_id}-region-{sanitize_id(gene)}",
            "meta": {
                "profile": [
                    "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/region-studied"
                ],
                "tag": [{"system": "http://terminology.kemkes.go.id/sp",
                         "code": "genomics", "display": "Genomics"}]
            },
            "text": {"status": "generated", "div": div_text},
            "status": "final",
            "category": [
                {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/observation-category",
                             "code": "laboratory", "display": "Laboratory"}]},
                {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                             "code": "GE", "display": "Genetics"}]}
            ],
            "code": {
                "coding": [{"system": "http://loinc.org", "code": "53041-0",
                            "display": "DNA region of interest panel"}],
                "text": f"Region studied: {gene}"
            },
            "subject": {"reference": f"Patient/{sample_id}-patient"},
            "specimen": {"reference": f"Specimen/{sample_id}-specimen"},
            "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
            "performer": [{"reference": f"Organization/{org_id}"}],
            "component": components
        })

    return observations


parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help='Path to VCF file')
parser.add_argument('--output', required=True, help='Path to output FHIR JSON')
parser.add_argument('--lineage_dir', help='Directory containing lineage JSON files')
parser.add_argument('--coverage_dir', help='Directory containing <sample>.coverage.json files')
parser.add_argument('--organization_metadata', default='', help='Path to organization_metadata CSV/Excel file')
parser.add_argument('--gene_map', default='', help='Optional TSV of gene<TAB>ncbi_gene_id to extend the built-in map')
parser.add_argument('--pipeline_version', default='unknown', help='Pipeline version, recorded in the bundle')
parser.add_argument('--coverage_min_depth', type=int, default=10)
parser.add_argument('--coverage_min_breadth', type=float, default=0.95)
parser.add_argument('--legacy-annotation-table', dest='legacy_annotation_table',
                    action='store_true',
                    help='Tolerate a one-row-per-pair annotation table by applying the most '
                         'severe grade to every listed drug. Over-calls resistance; prefer '
                         'rebuilding the table with scripts/rebuild_annotation_table.sh.')
args = parser.parse_args()

org_data = {}
if args.organization_metadata and os.path.exists(args.organization_metadata):
    org_data = load_organization_metadata(args.organization_metadata)
org_id = org_data.get('org_id', '100007732') if org_data else '100007732'

default_gene_map = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'data', 'gene_ncbi_map.tsv')
load_external_gene_map(args.gene_map or default_gene_map)

lineage_data = {}
if args.lineage_dir:
    lineage_data = load_lineage_data(args.lineage_dir)

coverage_data = {}
if args.coverage_dir:
    coverage_data = load_coverage_data(args.coverage_dir)

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

sample_coverage_info = None
for variation in file_sample_id_variations:
    if variation in coverage_data:
        sample_coverage_info = coverage_data[variation]
        break

if sample_coverage_info and sample_coverage_info.get('coverage_assessed'):
    assessable = sum(1 for d in sample_coverage_info.get('drugs', {}).values()
                     if d.get('assessable'))
    print(f"Coverage data found: {assessable} drugs assessable")
else:
    reason = (sample_coverage_info or {}).get('reason', 'no coverage file for this sample')
    print(f"No usable coverage data ({reason}); "
          f"non-resistant drugs will be reported as INDETERMINATE, not susceptible")

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
    skipped_records = []
    legacy_collapse_count = 0

    resistant_drugs_detected = set()

    for idx, record in enumerate(variants):
        try:
            file_sample_id = extract_sample_id_from_filename(args.input)[0]
            sample_id = f"{file_sample_id}-patient"
            pos = str(record.POS)
            ref = record.REF
            alt = str(record.ALT[0])

            gene = get_variant_value(record, 'GENE')
            effect = get_variant_value(record, 'EFFECT')
            genome_position = get_variant_value(record, 'GENOME_POSITION')
            depth_val = get_variant_value(record, 'DP', None)

            try:
                drug_pairs = resolve_variant_drug_pairs(
                    record, allow_legacy=args.legacy_annotation_table)
                if args.legacy_annotation_table and len(drug_pairs) > 1:
                    grades_seen = get_variant_list(record, 'WHO_CLASSIFICATION')
                    if len(grades_seen) == 1:
                        legacy_collapse_count += 1
            except AnnotationError as e:
                print(f"\nFATAL: {e}\n")
                sys.exit(2)

            drugs = [d for d, _g, _v in drug_pairs]
            grades = [g for _d, g, _v in drug_pairs]
            variant_ids = [v for _d, _g, v in drug_pairs if v]

            for drug, grade, _variant_id in drug_pairs:
                if is_resistance_grade(grade):
                    d_lower = drug.lower().strip()
                    resistant_drugs_detected.add(d_lower)
                    if 'fluoroquinolone' in d_lower:
                        resistant_drugs_detected.add('moxifloxacin')
                        resistant_drugs_detected.add('levofloxacin')
                        resistant_drugs_detected.add('ciprofloxacin')
                        resistant_drugs_detected.add('ofloxacin')

            if gene != "unknown" or drugs or grades:
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
                    "code": {"coding": [{"system": "http://loinc.org", "code": "48013-7",
                                         "display": "Genomic reference sequence ID"}]},
                    "valueCodeableConcept": {
                        "coding": [{"system": "http://www.ncbi.nlm.nih.gov/refseq",
                                    "code": "NC_000962.3", "display": "NC_000962.3"}],
                        "text": "NC_000962.3"
                    }
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

            allele_fraction = compute_allele_fraction(record)
            if allele_fraction is not None:
                components.append({
                    "code": {"coding": [{"system": "http://loinc.org", "code": "81258-6",
                                         "display": "Sample variant allelic frequency [NFr]"}]},
                    "valueQuantity": {
                        "value": round(allele_fraction, 4),
                        "system": "http://unitsofmeasure.org",
                        "code": "1"
                    }
                })

            if record.QUAL is not None:
                components.append({
                    "code": {
                        "coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                    "code": "variant-quality", "display": "Variant call quality (QUAL)"}],
                        "text": "Variant call quality (QUAL)"
                    },
                    "valueQuantity": {
                        "value": round(record.QUAL, 2),
                        "system": "http://unitsofmeasure.org",
                        "code": "1"
                    }
                })

            mapping_quality = get_variant_value(record, 'MQ', default=None)
            if mapping_quality is not None and mapping_quality != 'unknown':
                try:
                    components.append({
                        "code": {
                            "coding": [{"system": "http://hl7.org/fhir/uv/genomics-reporting/CodeSystem/tbd-codes-cs",
                                        "code": "mapping-quality", "display": "Mapping quality (MQ)"}],
                            "text": "Mapping quality (MQ)"
                        },
                        "valueQuantity": {
                            "value": round(float(mapping_quality), 2),
                            "system": "http://unitsofmeasure.org",
                            "code": "1"
                        }
                    })
                except ValueError:
                    pass

            if gene != "unknown":
                components.append({
                    "code": {"coding": [{"system": "http://loinc.org", "code": "48018-6", "display": "Gene studied [ID]"}]},
                    "valueCodeableConcept": build_gene_codeable_concept(gene)
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

            for drug, grade, _variant_id in drug_pairs:
                drug_snomed = get_drug_snomed_mapping(drug)
                if drug_snomed:
                    drug_value = {
                        "coding": [{
                            "system": "http://snomed.info/sct",
                            "code": drug_snomed['code'],
                            "display": drug_snomed['display']
                        }],
                        "text": drug
                    }
                else:
                    drug_value = {"text": drug}

                components.append({
                    "code": {"coding": [{"system": "http://loinc.org", "code": "51963-7",
                                         "display": "Medication assessed [Identifier]"}]},
                    "valueCodeableConcept": drug_value
                })

                who_coding = get_who_classification_coding(grade)
                if who_coding:
                    significance_value = {
                        "coding": [{
                            "system": who_coding['system'],
                            "code": who_coding['code'],
                            "display": who_coding['display']
                        }],
                        "text": f"{grade} ({drug})"
                    }
                else:
                    significance_value = {"text": f"{grade} ({drug})"}

                components.append({
                    "code": {"coding": [{"system": "http://loinc.org", "code": "53037-8",
                                         "display": "Genetic variation clinical significance [Imp]"}]},
                    "valueCodeableConcept": significance_value
                })

            for variant_id in dict.fromkeys(variant_ids):
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
                        }],
                        "text": variant_id
                    }
                })

            try:
                position_value = int(pos)
                components.append({
                    "code": {"coding": [{"system": "http://loinc.org", "code": "81254-5", "display": "Variant exact start-end"}]},
                    "valueRange": {
                        "low": {
                            "value": position_value
                        },
                        "high": {
                            "value": position_value + max(len(ref), 1) - 1
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
            if variant_ids:
                annotation_text += f" - {', '.join(dict.fromkeys(variant_ids))}"
            if drug_pairs:
                per_drug = "; ".join(f"{d}: {g}" for d, g, _v in drug_pairs)
                annotation_text += f" - WHO classification per drug - {per_drug}"

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
                "performer": [{"reference": f"Organization/{org_id}"}],
                "component": components
            }

            bundles[file_sample_id].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
                "resource": observation
            })

            variant_count += 1
            if variant_count % 100 == 0:
                print(f"Processed {variant_count} variants")

        except SystemExit:
            raise
        except Exception as e:
            skipped_records.append((idx, f"{record.CHROM}:{record.POS}", str(e)))
            print(f"Error processing variant {idx} at {record.CHROM}:{record.POS}: {e}")
            continue

    primary_sample_id = extract_sample_id_from_filename(args.input)[0]

    panel_drug_keys = {drug_key for drug_key, _c, _d in get_drug_panel_config()}
    coverage_assessed = bool(sample_coverage_info and sample_coverage_info.get('coverage_assessed'))
    coverage_drugs = (sample_coverage_info or {}).get('drugs', {})
    coverage_drugs_lower = {k.lower(): v for k, v in coverage_drugs.items()}

    panel_components = []
    panel_config = get_drug_panel_config()
    panel_state_counts = {'resistant': 0, 'no_resistance_detected': 0, 'indeterminate': 0}

    for drug_key, code, display in panel_config:
        is_resistant = drug_key in resistant_drugs_detected

        if is_resistant:
            value_code, value_display = "LA6676-6", "Resistant"
            state = 'resistant'
        else:
            drug_coverage = coverage_drugs_lower.get(drug_key)
            assessable = bool(coverage_assessed and drug_coverage
                              and drug_coverage.get('assessable'))
            if assessable:
                value_code, value_display = "LA24225-7", "Susceptible"
                state = 'no_resistance_detected'
            else:
                value_code, value_display = "LA9663-1", "Indeterminate"
                state = 'indeterminate'

        panel_state_counts[state] += 1

        component = {
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
        }

        if state == 'indeterminate':
            drug_coverage = coverage_drugs_lower.get(drug_key) or {}
            if not coverage_assessed:
                reason = (sample_coverage_info or {}).get(
                    'reason', 'no coverage data available for this sample')
            elif drug_coverage.get('inadequate_loci'):
                reason = ("insufficient coverage at: "
                          + ", ".join(drug_coverage['inadequate_loci']))
            elif drug_coverage.get('loci_without_coverage_data'):
                reason = ("no coverage data for: "
                          + ", ".join(drug_coverage['loci_without_coverage_data']))
            else:
                reason = 'resistance loci for this drug are not represented in the coverage targets'
            component["valueCodeableConcept"]["text"] = (
                f"Indeterminate - no resistance-associated mutation detected, but "
                f"susceptibility cannot be confirmed ({reason})"
            )

        panel_components.append(component)

    panel_summary = (f"{panel_state_counts['resistant']} resistant, "
                     f"{panel_state_counts['no_resistance_detected']} no resistance detected, "
                     f"{panel_state_counts['indeterminate']} indeterminate")
    print(f"Susceptibility panel: {panel_summary}")

    panel_observation = {
        "resourceType": "Observation",
        "id": f"{primary_sample_id}-susceptibility-panel",
        "meta": {
            "profile": ["http://hl7.org/fhir/StructureDefinition/Observation"],
        },
        "text": {
            "status": "generated",
            "div": (f"<div xmlns=\"http://www.w3.org/1999/xhtml\">Mycobacterial susceptibility "
                    f"panel for {primary_sample_id}: {panel_summary}. Indeterminate means no "
                    f"resistance-associated mutation was detected but coverage of the relevant "
                    f"loci could not be confirmed; it is not a susceptible result.</div>")
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
        "performer": [{"reference": f"Organization/{org_id}"}],
        "component": panel_components
    }

    bundles[primary_sample_id].append({
        "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
        "resource": panel_observation
    })

    for region_obs in build_region_studied_observations(
            primary_sample_id, sample_coverage_info, org_id, panel_drug_keys):
        bundles[primary_sample_id].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4()).lower()}",
            "resource": region_obs
        })

    if sample_lineage_info:
        lineage = sample_lineage_info.get('lineage', 'unknown')
        family = sample_lineage_info.get('family', 'Unknown')

        if lineage != 'unknown':
            lineage_obs_id = f"{primary_sample_id}-lineage"

            confidence = sample_lineage_info.get('confidence')
            score = sample_lineage_info.get('score')
            matched_snps = sample_lineage_info.get('matched_snps')
            total_snps = sample_lineage_info.get('total_snps')

            lineage_text = f"Mycobacterial Lineage: {lineage}"
            if family != 'Unknown':
                lineage_text += f" ({family})"
            if matched_snps is not None and total_snps:
                lineage_text += f" - {matched_snps}/{total_snps} barcode SNPs"
            if confidence:
                lineage_text += f", confidence {confidence}"

            div_text = f"<div xmlns=\"http://www.w3.org/1999/xhtml\">{lineage_text}</div>"

            lineage_components = []

            if family and family != 'Unknown':
                lineage_components.append({
                    "code": {"coding": [{"system": "http://terminology.spheres.id/CodeSystem/mtb-lineage-attribute",
                                         "code": "lineage-family", "display": "Lineage family"}],
                             "text": "Lineage family"},
                    "valueString": family
                })

            if confidence:
                lineage_components.append({
                    "code": {"coding": [{"system": "http://terminology.spheres.id/CodeSystem/mtb-lineage-attribute",
                                         "code": "lineage-confidence", "display": "Lineage call confidence"}],
                             "text": "Lineage call confidence"},
                    "valueString": str(confidence)
                })

            if score is not None:
                try:
                    lineage_components.append({
                        "code": {"coding": [{"system": "http://terminology.spheres.id/CodeSystem/mtb-lineage-attribute",
                                             "code": "lineage-score", "display": "Lineage barcode score"}],
                                 "text": "Lineage barcode score"},
                        "valueQuantity": {"value": float(score),
                                          "system": "http://unitsofmeasure.org", "code": "1"}
                    })
                except (TypeError, ValueError):
                    pass

            if matched_snps is not None:
                try:
                    lineage_components.append({
                        "code": {"coding": [{"system": "http://terminology.spheres.id/CodeSystem/mtb-lineage-attribute",
                                             "code": "barcode-snps-matched", "display": "Barcode SNPs matched"}],
                                 "text": "Barcode SNPs matched"},
                        "valueQuantity": {"value": int(matched_snps),
                                          "system": "http://unitsofmeasure.org", "code": "1"}
                    })
                except (TypeError, ValueError):
                    pass

            if total_snps is not None:
                try:
                    lineage_components.append({
                        "code": {"coding": [{"system": "http://terminology.spheres.id/CodeSystem/mtb-lineage-attribute",
                                             "code": "barcode-snps-total", "display": "Barcode SNPs examined"}],
                                 "text": "Barcode SNPs examined"},
                        "valueQuantity": {"value": int(total_snps),
                                          "system": "http://unitsofmeasure.org", "code": "1"}
                    })
                except (TypeError, ValueError):
                    pass

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
                    "coding": [
                        {
                            "system": "http://terminology.spheres.id/CodeSystem/mtb-lineage",
                            "code": lineage,
                            "display": f"M. tuberculosis {lineage}"
                        },
                        {
                            "system": "http://tb-lineage.org",
                            "code": lineage,
                            "display": f"TB Lineage {lineage}"
                        }
                    ],
                    "text": f"Lineage {lineage}"
                },
                "subject": {"reference": f"Patient/{primary_sample_id}-patient"},
                "specimen": {"reference": f"Specimen/{primary_sample_id}-specimen"},
                "effectiveDateTime": datetime.now(timezone.utc).isoformat(),
                "performer": [{"reference": f"Organization/{org_id}"}]
            }

            if lineage_components:
                lineage_observation["component"] = lineage_components

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

    print(f"\nSummary for {primary_sample_id}:")
    print(f"  variants written        : {variant_count}")
    print(f"  annotated variants      : {successful_annotations}")
    print(f"  records skipped (error) : {len(skipped_records)}")
    if skipped_records:
        for idx, location, message in skipped_records[:10]:
            print(f"      - record {idx} at {location}: {message}")
        if len(skipped_records) > 10:
            print(f"      ... and {len(skipped_records) - 10} more")
    print(f"  resistant drugs detected: "
          f"{', '.join(sorted(resistant_drugs_detected)) if resistant_drugs_detected else 'none'}")
    print(f"  susceptibility panel    : {panel_summary}")

    if legacy_collapse_count:
        print(f"  WARNING: {legacy_collapse_count} multi-drug alleles had a single "
              f"collapsed grade applied to every listed drug (--legacy-annotation-table). "
              f"Rebuild the annotation table to remove this ambiguity.")

    if _unmapped_genes:
        print(f"  genes without an NCBI GeneID: {len(_unmapped_genes)} "
              f"(emitted as text only, no fabricated coding)")
        print(f"      {', '.join(sorted(_unmapped_genes))}")
        if not _external_gene_map_loaded:
            print(f"      Populate data/gene_ncbi_map.tsv (gene<TAB>ncbi_gene_id) to code these.")

    try:
        os.remove(fixed_vcf)
    except:
        pass

except SystemExit:
    raise
except Exception as e:
    import traceback
    print(f"Error: {e}")
    traceback.print_exc()
    sys.exit(1)
