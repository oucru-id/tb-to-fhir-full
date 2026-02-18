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
                    
                    variant_data = {
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'gene': info.get('GENE', 'unknown'),
                        'drug': info.get('DRUG', 'unknown'),
                        'effect': info.get('EFFECT', 'unknown'),
                        'who_classification': info.get('WHO_CLASSIFICATION', 'unknown'),
                        'variant_id': info.get('VARIANT_ID', 'unknown'),
                        'genome_position': info.get('GENOME_POSITION', 'unknown')
                    }
                    
                    variants.append(variant_data)
                        
    except Exception as e:
        debug_print(f"Error parsing VCF file: {e}")
    
    return variants

def generate_summary_report(sample_id, variants, lineage_info, output_dir):
    output_file = os.path.join(output_dir, f"{sample_id}.summary_report.txt")
    
    resistance_variants = []
    interim_variants = []
    not_assoc_variants = []
    
    for variant in variants:
        who_class = variant['who_classification']
        if who_class != 'unknown':        
            if 'Assoc w R' in who_class and 'Interim' not in who_class:
                resistance_variants.append(variant)
            elif 'Assoc w R - Interim' in who_class:
                interim_variants.append(variant)
            else:
                not_assoc_variants.append(variant)
    
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
        f.write(f"Confirmed resistance: {len(resistance_variants)}\n")
        f.write(f"Interim resistance: {len(interim_variants)}\n")
        f.write(f"Not assoc w R: {len(not_assoc_variants)}\n\n")
        

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
        
        f.write(f"\nCLINICAL SUMMARY:\n")
        f.write("=" * 30 + "\n")
        
        if lineage_info:
            lineage = lineage_info.get('lineage', 'Unknown')
            family = lineage_info.get('family', 'Unknown')
            f.write(f"Lineage: {lineage} ({family})\n")
        
        if resistance_variants:
            f.write(f"Drug resistance: DETECTED ({len(resistance_variants)} variants)\n")
            resistant_drugs = set(v['drug'] for v in resistance_variants if v['drug'] != 'unknown')
            if resistant_drugs:
                f.write(f"Resistant to: {', '.join(sorted(resistant_drugs))}\n")
            f.write(f"Clinical action: REQUIRED\n")
        else:
            f.write(f"Drug resistance: NOT DETECTED\n")
            f.write(f"Clinical action: Standard treatment\n")
        
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