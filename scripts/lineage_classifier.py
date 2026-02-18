#!/usr/bin/env python3

import json
import argparse
import sys
from collections import defaultdict
from datetime import datetime


def load_barcode_bed(barcode_bed_path):
    lineage_snps = defaultdict(list)
    lineage_info = {}
    
    try:
        with open(barcode_bed_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    lineage = parts[3]
                    ref_alt = parts[4]  
                    lineage_family = parts[5] if len(parts) > 5 else "Unknown"
                    
                    pos = end
                    
                    snp_key = f"{chrom}:{pos}:{ref_alt}"
                    lineage_snps[lineage].append({
                        'chrom': chrom,
                        'pos': pos,
                        'alt': ref_alt,
                        'snp_key': snp_key
                    })
                    
                    if lineage not in lineage_info:
                        lineage_info[lineage] = {
                            'family': lineage_family,
                            'snp_count': 0
                        }
                    lineage_info[lineage]['snp_count'] += 1
        
        return lineage_snps, lineage_info
        
    except Exception as e:
        return {}, {}

def extract_sample_snps(vcf_path):
    sample_snps = {}
    variant_count = 0
    
    try:
        with open(vcf_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                variant_count += 1
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    chrom = parts[0]
                    pos = int(parts[1])
                    ref = parts[3]
                    alt = parts[4]
                    
                    alts = alt.split(',')
                    for alt_allele in alts:
                        snp_key = f"{chrom}:{pos}:{alt_allele}"
                        sample_snps[snp_key] = {
                            'chrom': chrom,
                            'pos': pos,
                            'ref': ref,
                            'alt': alt_allele
                        }
        
        return sample_snps
        
    except Exception as e:
        return {}

def classify_lineage(sample_snps, lineage_snps, lineage_info):
    lineage_scores = {}
    
    for lineage_name, snps in lineage_snps.items():
        total_snps = len(snps)
        matched_snps = 0
        matched_details = []
        
        for snp in snps:
            snp_key = snp['snp_key']
            if snp_key in sample_snps:
                matched_snps += 1
                matched_details.append({
                    'pos': snp['pos'],
                    'alt': snp['alt'],
                    'snp_key': snp_key
                })
        
        if total_snps > 0:
            score = matched_snps / total_snps
            lineage_scores[lineage_name] = {
                'score': score,
                'matched': matched_snps,
                'total': total_snps,
                'percentage': round(score * 100, 2),
                'family': lineage_info.get(lineage_name, {}).get('family', 'Unknown'),
                'matched_details': matched_details
            }
    
    if lineage_scores:
        best_lineage = max(lineage_scores.keys(), 
                          key=lambda x: (lineage_scores[x]['score'], lineage_scores[x]['matched']))
        best_score = lineage_scores[best_lineage]
        
        confidence = 'low'
        if best_score['score'] >= 0.8 and best_score['matched'] >= 3:
            confidence = 'high'
        elif best_score['score'] >= 0.6 and best_score['matched'] >= 2:
            confidence = 'medium'
        
        return {
            'lineage': best_lineage,
            'confidence': confidence,
            'score': best_score['score'],
            'matched_snps': best_score['matched'],
            'total_snps': best_score['total'],
            'percentage': best_score['percentage'],
            'family': best_score['family'],
            'matched_details': best_score['matched_details'],
            'all_scores': {k: v for k, v in lineage_scores.items() if v['matched'] > 0}
        }
    
    return {
        'lineage': 'unknown',
        'confidence': 'low',
        'score': 0,
        'matched_snps': 0,
        'total_snps': 0,
        'percentage': 0,
        'family': 'Unknown',
        'matched_details': [],
        'all_scores': {}
    }

def main():
    parser = argparse.ArgumentParser(description='Classify TB lineage from barcode SNPs using BED format')
    parser.add_argument('--vcf', required=True, help='VCF file with barcode variants')
    parser.add_argument('--barcode', required=True, help='Barcode BED file')
    parser.add_argument('--sample_id', required=True, help='Sample ID')
    parser.add_argument('--output', required=True, help='Output JSON file')
    
    args = parser.parse_args()
    
    lineage_snps, lineage_info = load_barcode_bed(args.barcode)
    if not lineage_snps:
        result = {
            'lineage': 'unknown',
            'confidence': 'low',
            'score': 0,
            'matched_snps': 0,
            'total_snps': 0,
            'percentage': 0,
            'family': 'Unknown',
            'sample_id': args.sample_id,
            'analysis_date': datetime.now().isoformat(),
            'error': 'Failed to load lineage barcode data'
        }
    else:
        sample_snps = extract_sample_snps(args.vcf)
        
        result = classify_lineage(sample_snps, lineage_snps, lineage_info)
        
        result['sample_id'] = args.sample_id
    
    with open(args.output, 'w') as f:
        json.dump(result, f, indent=2)

if __name__ == '__main__':
    main()