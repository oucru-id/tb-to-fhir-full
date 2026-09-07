#!/usr/bin/env python3

import argparse
import os
import sys
import json
import glob
from datetime import datetime
import gzip
import re

def debug_print(message):
    print(f"DEBUG: {message}", file=sys.stderr)

def extract_sample_id(filename):
    basename = os.path.basename(filename)
    
    basename = basename.replace('.annotated_variants.vcf.gz', '')
    basename = basename.replace('.annotated_variants.vcf', '')
    basename = basename.replace('.vcf.gz', '')
    basename = basename.replace('.vcf', '')
    
    basename = basename.replace('_ont.fastq', '')
    basename = basename.replace('_illumina', '')
    
    core_id = re.sub(r'_ont\.fastq$', '', basename)
    core_id = re.sub(r'_illumina$', '', core_id)
    
    return core_id

def _split_info(value):
    """Split a Number=. INFO value into its elements."""
    if not value or value is True or str(value).lower() == 'unknown':
        return []
    parts = [v.strip() for v in str(value).split(',')]
    return [p for p in parts if p and p.lower() != 'unknown']


def parse_vcf_info(info_str):
    info_dict = {}
    if info_str != '.':
        for item in info_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info_dict[key] = value
            else:
                info_dict[item] = True
    return info_dict

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
                
                clean_sample_id = sample_id.replace('_ont.fastq', '').replace('_illumina', '')
                
                lineage_data[sample_id] = lineage_info
                lineage_data[clean_sample_id] = lineage_info
                                
        except Exception as e:
            debug_print(f"Error loading lineage file {lineage_file}: {e}")
    
    return lineage_data

def parse_annotated_vcf(vcf_file):
    variants = []
    
    
    if vcf_file.endswith('.gz'):
        file_handle = gzip.open(vcf_file, 'rt')
    else:
        file_handle = open(vcf_file, 'r')
    
    try:
        with file_handle as f:
            line_count = 0
            for line in f:
                line_count += 1
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    info = parse_vcf_info(fields[7])
                    
                    if len(variants) < 5:
                        debug_print(f"Variant {len(variants)+1}: pos={pos}, ref={ref}, alt={alt}")
                        debug_print(f"  INFO keys: {list(info.keys())}")
                        debug_print(f"  WHO_CLASSIFICATION: {info.get('WHO_CLASSIFICATION', 'NOT FOUND')}")
                        debug_print(f"  GENE: {info.get('GENE', 'NOT FOUND')}")
                        debug_print(f"  DRUG: {info.get('DRUG', 'NOT FOUND')}")
                    
                    drugs = _split_info(info.get('DRUG'))
                    grades = _split_info(info.get('WHO_CLASSIFICATION'))
                    variant_ids = _split_info(info.get('VARIANT_ID'))

                    drug_pairs = []
                    if drugs and len(grades) == len(drugs):
                        for i, drug in enumerate(drugs):
                            drug_pairs.append({
                                'drug': drug,
                                'who_classification': grades[i],
                                'variant_id': variant_ids[i] if i < len(variant_ids) else (
                                    variant_ids[0] if len(variant_ids) == 1 else 'unknown'),
                            })
                    elif drugs:
                        debug_print(f"WARNING {chrom}:{pos} {ref}>{alt}: {len(drugs)} drugs but "
                                    f"{len(grades)} grades -- annotation table needs rebuilding "
                                    f"(scripts/rebuild_annotation_table.sh)")
                        for i, drug in enumerate(drugs):
                            drug_pairs.append({
                                'drug': drug,
                                'who_classification': grades[0] if grades else 'unknown',
                                'variant_id': variant_ids[0] if variant_ids else 'unknown',
                                'ambiguous': True,
                            })

                    variant_data = {
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'gene': info.get('GENE', 'unknown'),
                        'effect': info.get('EFFECT', 'unknown'),
                        'genome_position': info.get('GENOME_POSITION', 'unknown'),
                        'drug_pairs': drug_pairs,
                    }

                    variants.append(variant_data)
                        
    except Exception as e:
        debug_print(f"Error parsing VCF file: {e}")
    
    return variants

def generate_summary_report(sample_id, variants, lineage_info, output_dir):
    output_file = os.path.join(output_dir, f"{sample_id}.summary_report.txt")
    
    resistance_variants = []
    interim_variants = []
    uncertain_variants = []
    not_assoc_variants = []
    ambiguous_pairs = 0
    graded_pair_count = 0

    for variant in variants:
        for pair in variant.get('drug_pairs', []):
            who_class = pair.get('who_classification', 'unknown')
            if who_class == 'unknown':
                continue

            graded_pair_count += 1
            if pair.get('ambiguous'):
                ambiguous_pairs += 1

            finding = dict(variant)
            finding.update(pair)

            normalized = who_class.strip().lower()
            if normalized == 'assoc w r':
                resistance_variants.append(finding)
            elif normalized == 'assoc w r - interim':
                interim_variants.append(finding)
            elif normalized == 'uncertain significance':
                uncertain_variants.append(finding)
            else:
                not_assoc_variants.append(finding)

    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write(f"TB GENOMIC ANALYSIS SUMMARY REPORT - {sample_id}\n")
        f.write("="*80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("LINEAGE CLASSIFICATION:\n")
        f.write("=" * 50 + "\n")
        if lineage_info:
            lineage = lineage_info.get('lineage', 'Unknown')
            family = lineage_info.get('family', 'Unknown')
            
            f.write(f"Lineage: {lineage}\n")
            f.write(f"Family: {family}\n")
        else:
            f.write("LINEAGE: Not determined\n")
        
        f.write(f"\nRESISTANCE ANALYSIS:\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total variants analyzed: {len(variants)}\n")
        f.write(f"Graded variant-drug pairs: {graded_pair_count}\n")
        f.write(f"Confirmed resistance (WHO group 1): {len(resistance_variants)}\n")
        f.write(f"Interim resistance (WHO group 2): {len(interim_variants)}\n")
        f.write(f"Uncertain significance (WHO group 3): {len(uncertain_variants)}\n")
        f.write(f"Not assoc w R (WHO groups 4-5): {len(not_assoc_variants)}\n\n")

        if ambiguous_pairs:
            f.write("WARNING: the annotation table has one row per (variant, drug) pair, so\n")
            f.write(f"         {ambiguous_pairs} pair(s) share a single collapsed WHO grade and\n")
            f.write("         their per-drug association is not reliable. Rebuild the table\n")
            f.write("         with scripts/rebuild_annotation_table.sh.\n\n")


        if resistance_variants:
            f.write("CONFIRMED RESISTANCE VARIANTS (WHO: Assoc w R):\n")
            f.write("=" * 60 + "\n")
            

            drug_groups = {}
            for variant in resistance_variants:
                drug = variant['drug']
                if drug not in drug_groups:
                    drug_groups[drug] = []
                drug_groups[drug].append(variant)
            
            for drug, drug_variants in sorted(drug_groups.items()):
                f.write(f"\nDRUG: {drug.upper()}\n")
                f.write("-" * 40 + "\n")
                

                for i, variant in enumerate(drug_variants, 1):
                    f.write(f"Variant {i}:\n")
                    f.write(f"  Position: {variant['pos']}\n")
                    f.write(f"  Gene: {variant['gene']}\n")
                    f.write(f"  Effect: {variant['effect']}\n")
                    if variant['variant_id'] != 'unknown':
                        f.write(f"  Variant ID: {variant['variant_id']}\n")
                    f.write(f"  Mutation: {variant['ref']} -> {variant['alt']}\n")
                    if i < len(drug_variants):  
                        f.write("  " + "-" * 30 + "\n")
        

        if interim_variants:
            f.write(f"\nINTERIM RESISTANCE VARIANTS (WHO: Assoc w R - Interim):\n")
            f.write("-" * 60 + "\n")
            f.write("These variants may be associated with resistance\n")
            
            drug_groups = {}
            for variant in interim_variants:
                drug = variant['drug']
                if drug not in drug_groups:
                    drug_groups[drug] = []
                drug_groups[drug].append(variant)
            
            for drug, drug_variants in sorted(drug_groups.items()):
                f.write(f"\nDRUG: {drug.upper()}\n")
                f.write("-" * 40 + "\n")
                
                for i, variant in enumerate(drug_variants, 1):
                    f.write(f"Variant {i}:\n")
                    f.write(f"  Position: {variant['pos']}\n")
                    f.write(f"  Gene: {variant['gene']}\n")
                    f.write(f"  Effect: {variant['effect']}\n")
                    if variant['variant_id'] != 'unknown':
                        f.write(f"  Variant ID: {variant['variant_id']}\n")
                    f.write(f"  Mutation: {variant['ref']} -> {variant['alt']}\n")
                    if i < len(drug_variants):
                        f.write("  " + "-" * 30 + "\n")
        
        if uncertain_variants:
            f.write(f"\nVARIANTS OF UNCERTAIN SIGNIFICANCE (WHO group 3):\n")
            f.write("-" * 60 + "\n")
            f.write("WHO group 3 is NOT a susceptible result. These variants have insufficient\n")
            f.write("or conflicting evidence and may still be associated with resistance.\n")

            drug_groups = {}
            for variant in uncertain_variants:
                drug_groups.setdefault(variant['drug'], []).append(variant)

            for drug, drug_variants in sorted(drug_groups.items()):
                f.write(f"\nDRUG: {drug.upper()}\n")
                f.write("-" * 40 + "\n")
                for i, variant in enumerate(drug_variants, 1):
                    f.write(f"Variant {i}: {variant['gene']} "
                            f"{variant.get('variant_id', 'unknown')} "
                            f"at {variant['pos']} ({variant['ref']} -> {variant['alt']})\n")

        f.write(f"\nCLINICAL SUMMARY:\n")
        f.write("=" * 30 + "\n")

        if lineage_info:
            lineage = lineage_info.get('lineage', 'Unknown')
            family = lineage_info.get('family', 'Unknown')
            f.write(f"Lineage: {lineage} ({family})\n")

        if resistance_variants or interim_variants:
            total = len(resistance_variants) + len(interim_variants)
            f.write(f"Drug resistance: DETECTED ({total} variant-drug findings)\n")
            resistant_drugs = set(v['drug'] for v in (resistance_variants + interim_variants)
                                  if v['drug'] != 'unknown')
            if resistant_drugs:
                f.write(f"Resistant to: {', '.join(sorted(resistant_drugs))}\n")
            f.write(f"Clinical action: REQUIRED\n")
        else:
            f.write(f"Drug resistance: NO RESISTANCE-ASSOCIATED MUTATION DETECTED\n")
            f.write(f"Note: this is not a susceptibility result. Susceptibility can only be\n")
            f.write(f"      inferred where coverage of the relevant resistance loci was\n")
            f.write(f"      confirmed; see the susceptibility panel for per-drug assessability.\n")

        if uncertain_variants:
            uncertain_drugs = sorted(set(v['drug'] for v in uncertain_variants
                                         if v['drug'] != 'unknown'))
            f.write(f"Uncertain-significance findings for: {', '.join(uncertain_drugs)}\n")

        genes_with_resistance = set(v['gene'] for v in resistance_variants if v['gene'] != 'unknown')
        if genes_with_resistance:
            f.write(f"Genes with resistance: {', '.join(sorted(genes_with_resistance))}\n")

        f.write(f"\n" + "="*80 + "\n")
        f.write("END OF SUMMARY REPORT\n")
        f.write("="*80 + "\n")

def main():
    parser = argparse.ArgumentParser(description='Generate TB sample summary reports')
    parser.add_argument('--annotated_vcf', required=True, help='Annotated VCF file')
    parser.add_argument('--lineage_dir', help='Directory containing lineage JSON files')
    parser.add_argument('--output_dir', default='.', help='Output directory')
    args = parser.parse_args()
    
    sample_id = extract_sample_id(args.annotated_vcf)
    
    variants = parse_annotated_vcf(args.annotated_vcf)
    
    lineage_data = load_lineage_data(args.lineage_dir) if args.lineage_dir else {}
    
    lineage_info = None
    for potential_id in [sample_id, sample_id.replace('_ont.fastq', ''), sample_id.replace('_illumina', '')]:
        if potential_id in lineage_data:
            lineage_info = lineage_data[potential_id]
            break
    
    if not lineage_info:
        debug_print(f"No lineage data found for sample {sample_id}")
    
    generate_summary_report(sample_id, variants, lineage_info, args.output_dir)
    
if __name__ == '__main__':
    main()
