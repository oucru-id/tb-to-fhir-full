#!/usr/bin/env python3

import pandas as pd
import argparse
import gzip
import os
import re
import sys
from collections import OrderedDict

def cross_check_coordinates(genomic_df, vcf_path):

    from collections import defaultdict

    vcf_alleles = defaultdict(set)
    records = 0
    opener = gzip.open if str(vcf_path).endswith('.gz') else open
    with opener(vcf_path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                continue
            records += 1
            info = fields[7]
            if 'graded_variant=' not in info:
                continue
            graded = info.split('graded_variant=', 1)[1].split(';')[0]
            key = (fields[0], int(fields[1]), fields[3], fields[4])
            vcf_alleles[key].update(g.strip() for g in graded.split('&') if g.strip())

    sheet_alleles = defaultdict(set)
    for row in genomic_df.itertuples(index=False):
        key = (str(row.CHROM), int(row.POS), str(row.REF), str(row.ALT))
        sheet_alleles[key].add(str(row.variant).strip())

    only_vcf = set(vcf_alleles) - set(sheet_alleles)
    only_sheet = set(sheet_alleles) - set(vcf_alleles)
    differing = [k for k in (set(vcf_alleles) & set(sheet_alleles))
                 if vcf_alleles[k] != sheet_alleles[k]]

    print(f"Cross-check against {os.path.basename(str(vcf_path))}:")
    print(f"  VCF records / unique alleles : {records} / {len(vcf_alleles)}")
    print(f"  spreadsheet unique alleles   : {len(sheet_alleles)}")
    print(f"  alleles only in VCF          : {len(only_vcf)}")
    print(f"  alleles only in spreadsheet  : {len(only_sheet)}")
    print(f"  differing graded-variant sets: {len(differing)}")

    if only_vcf or only_sheet or differing:
        print("  MISMATCH: the spreadsheet and the coordinate VCF disagree. They are "
              "probably from different catalogue releases; resolve before using either.")
        for key in list(only_vcf)[:3]:
            print(f"    only in VCF        : {key} -> {sorted(vcf_alleles[key])}")
        for key in list(only_sheet)[:3]:
            print(f"    only in spreadsheet: {key} -> {sorted(sheet_alleles[key])}")
        for key in differing[:3]:
            print(f"    differs            : {key}")
            print(f"        VCF   : {sorted(vcf_alleles[key])}")
            print(f"        sheet : {sorted(sheet_alleles[key])}")
        return False

    print("  OK: coordinate sources agree exactly.")
    return True


def main():
    parser = argparse.ArgumentParser(description="Convert WHO-UCN-TB-2023 Excel file to TSV")
    parser.add_argument('--input', required=True, help='WHO-UCN-TB-2023.xlsx file')
    parser.add_argument('--output', default='enhanced_annotation_table.tsv', help='Output TSV file')
    parser.add_argument('--header', default='enhanced_annotations_header.txt', help='Output header file')
    parser.add_argument('--coordinates_vcf', default=None,
                        help="WHO's own coordinate VCF (Genomic_coordinates_*.vcf.gz). When "
                             "given, the extracted coordinates are verified against it and "
                             "any disagreement aborts the build.")
    args = parser.parse_args()

    try:
        
        genomic_df = pd.read_excel(args.input, sheet_name="Genomic_coordinates")
        print(f"Loaded {len(genomic_df)} genomic variants")
        
        catalogue_df = pd.read_excel(args.input, sheet_name="Catalogue_master_file", header=2)
        print(f"Loaded {len(catalogue_df)} graded variants")
        
        genomic_df = genomic_df.rename(columns={
            'chromosome': 'CHROM',
            'position': 'POS',
            'reference_nucleotide': 'REF',
            'alternative_nucleotide': 'ALT'
        })
        if 'CHROM' not in genomic_df.columns:
            genomic_df['CHROM'] = 'NC_000962.3'

        distinct_refs = sorted(set(genomic_df['CHROM'].dropna().astype(str)))
        print(f"Genomic reference(s) in coordinates sheet: {', '.join(distinct_refs)}")
        if distinct_refs != ['NC_000962.3']:
            print("Warning: coordinates are not exclusively NC_000962.3; the pipeline "
                  "reference must match or variants will not be matched.")

        if args.coordinates_vcf:
            if not os.path.exists(args.coordinates_vcf):
                print(f"ERROR: coordinate VCF not found: {args.coordinates_vcf}")
                sys.exit(1)
            if not cross_check_coordinates(genomic_df, args.coordinates_vcf):
                sys.exit(1)
        
        variant_col = 'variant'
        drug_col = 'drug'
        gene_col = 'gene'
        effect_col = 'effect'  
        
        who_classification_col = 'Initial confidence grading ALL dataset'
        
        if who_classification_col not in catalogue_df.columns:
            print(f"Warning: Column '{who_classification_col}' not found. Available columns:")
            for col in catalogue_df.columns:
                if 'confidence' in col.lower() or 'grading' in col.lower() or 'initial' in col.lower():
                    print(f"  - {col}")
            for alt_col in ['Initial confidence grading ALL dataset', 'confidence', 'grading', 'WHO_classification']:
                if alt_col in catalogue_df.columns:
                    who_classification_col = alt_col
                    break
        
        print(f"Using WHO classification column: {who_classification_col}")
        
        required_cols = [variant_col, drug_col, gene_col, effect_col, who_classification_col]
        
        existing_cols = []
        for col in required_cols:
            if col in catalogue_df.columns:
                existing_cols.append(col)
            else:
                print(f"Warning: Column '{col}' not found in catalogue")
        
        catalogue_subset = catalogue_df[existing_cols].copy()
        
        column_mapping = {
            variant_col: 'variant',
            drug_col: 'drug',
            gene_col: 'gene',
            effect_col: 'effect',
            who_classification_col: 'WHO_classification'
        }
        
        catalogue_subset = catalogue_subset.rename(columns=column_mapping)
        
        merged_df = pd.merge(
            genomic_df[['CHROM', 'POS', 'REF', 'ALT', 'variant']], 
            catalogue_subset, 
            on='variant', 
            how='left'
        )
        
        merged_df = merged_df.fillna('unknown')
        
        def standardize_effect(effect):
            if pd.isna(effect) or effect == 'unknown':
                return 'unknown'
            
            effect = str(effect).strip()
            
            effect_mapping = {
                'missense_variant': 'missense_variant',
                'frameshift': 'frameshift_variant',
                'inframe_deletion': 'inframe_deletion',
                'inframe_insertion': 'inframe_insertion',
                'feature_ablation': 'feature_ablation',
                'initiator_codon_variant': 'initiator_codon_variant',
                'LoF': 'loss_of_function_variant',
                'synonymous_variant': 'synonymous_variant',
                'stop_gained': 'stop_gained',
                'stop_lost': 'stop_lost',
                'splice_site_variant': 'splice_site_variant',
                'upstream_gene_variant': 'upstream_gene_variant',
                'downstream_gene_variant': 'downstream_gene_variant',
                'intergenic_variant': 'intergenic_variant',
                'intron_variant': 'intron_variant',
                'regulatory_region_variant': 'regulatory_region_variant'
            }
            
            return effect_mapping.get(effect, effect)
        
        merged_df['effect'] = merged_df['effect'].apply(standardize_effect)
        
        WHO_GRADE_BY_NUMBER = {
            '1': 'Assoc w R',
            '2': 'Assoc w R - Interim',
            '3': 'Uncertain significance',
            '4': 'Not assoc w R - Interim',
            '5': 'Not assoc w R',
        }
        WHO_GRADE_TEXTS = [
            'Not assoc w R - Interim',
            'Assoc w R - Interim',
            'Uncertain significance',
            'Not assoc w R',
            'Assoc w R',
        ]

        def standardize_who_classification(classification):
            if pd.isna(classification) or classification == 'unknown':
                return 'unknown'

            classification = str(classification).strip()

            leading = re.match(r'^\s*([1-5])\b', classification)
            if leading:
                return WHO_GRADE_BY_NUMBER[leading.group(1)]

            for text in WHO_GRADE_TEXTS:
                if text.lower() in classification.lower():
                    return text

            return classification

        merged_df['WHO_classification'] = merged_df['WHO_classification'].apply(standardize_who_classification)
        
        flat_df = pd.DataFrame({
            'CHROM': merged_df['CHROM'],
            'POS': merged_df['POS'],
            'REF': merged_df['REF'],
            'ALT': merged_df['ALT'],
            'GENE': merged_df.get('gene', 'unknown'),
            'DRUG': merged_df.get('drug', 'unknown'),
            'EFFECT': merged_df['effect'],
            'WHO_CLASSIFICATION': merged_df['WHO_classification'],
            'VARIANT_ID': merged_df['variant'],
            'GENOME_POSITION': merged_df['POS']
        })

        flat_df = flat_df.drop_duplicates().sort_values(['CHROM', 'POS', 'REF', 'ALT'])

        def _clean(value):
            return str(value).strip()

        def _is_known(value):
            return value and value.lower() != 'unknown'

        def _first_known(values):
            for value in values:
                text = _clean(value)
                if _is_known(text):
                    return text
            return 'unknown'

        print("Collapsing to one row per allele with parallel drug/grade arrays...")

        grouped = OrderedDict()
        for row in flat_df.itertuples(index=False):
            key = (row.CHROM, row.POS, row.REF, row.ALT)
            grouped.setdefault(key, []).append(row)

        collapsed_rows = []
        multi_drug_alleles = 0

        ungraded_pairs_dropped = 0

        for (chrom, pos, ref, alt), rows in grouped.items():
            graded = [r for r in rows
                      if _is_known(_clean(r.DRUG)) and _is_known(_clean(r.WHO_CLASSIFICATION))]
            ungraded_pairs_dropped += sum(
                1 for r in rows
                if _is_known(_clean(r.DRUG)) and not _is_known(_clean(r.WHO_CLASSIFICATION)))

            if not graded:
                collapsed_rows.append((
                    chrom, pos, ref, alt,
                    _first_known([r.GENE for r in rows]),
                    'unknown',
                    _first_known([r.EFFECT for r in rows]),
                    'unknown',
                    _first_known([r.VARIANT_ID for r in rows]),
                    rows[0].GENOME_POSITION,
                ))
                continue

            seen = set()
            pairs = []
            for r in graded:
                triple = (_clean(r.DRUG), _clean(r.WHO_CLASSIFICATION), _clean(r.VARIANT_ID))
                if triple not in seen:
                    seen.add(triple)
                    pairs.append(triple)
            pairs.sort(key=lambda t: (t[0], t[2]))

            for drug, grade, variant_id in pairs:
                if ',' in drug or ',' in grade or ',' in variant_id:
                    print(f"ERROR: comma inside a field at {chrom}:{pos} {ref}>{alt} "
                          f"would corrupt the parallel arrays: {drug!r} {grade!r} {variant_id!r}")
                    sys.exit(1)

            if len(pairs) > 1:
                multi_drug_alleles += 1

            collapsed_rows.append((
                chrom, pos, ref, alt,
                _first_known([r.GENE for r in graded]),
                ','.join(p[0] for p in pairs),
                _first_known([r.EFFECT for r in graded]),
                ','.join(p[1] for p in pairs),
                ','.join(p[2] for p in pairs),
                graded[0].GENOME_POSITION,
            ))

        final_df = pd.DataFrame(collapsed_rows, columns=[
            'CHROM', 'POS', 'REF', 'ALT', 'GENE', 'DRUG',
            'EFFECT', 'WHO_CLASSIFICATION', 'VARIANT_ID', 'GENOME_POSITION'
        ]).sort_values(['CHROM', 'POS', 'REF', 'ALT'])

        print(f"Alleles collapsed: {len(flat_df)} pair-rows -> {len(final_df)} allele-rows")
        print(f"  multi-drug alleles: {multi_drug_alleles}")
        if ungraded_pairs_dropped:
            print(f"  drug associations dropped for having no grade in "
                  f"'{who_classification_col}': {ungraded_pairs_dropped}")

        final_df.to_csv(args.output, sep='\t', index=False, header=False)
        print(f"Enhanced annotation table saved to {args.output}")
        print(f"Total variants: {len(final_df)}")

        header_content = """##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
##INFO=<ID=DRUG,Number=.,Type=String,Description="Associated drug(s)">
##INFO=<ID=EFFECT,Number=1,Type=String,Description="Predicted effect of mutation">
##INFO=<ID=WHO_CLASSIFICATION,Number=.,Type=String,Description="WHO resistance classification, one per DRUG, same order">
##INFO=<ID=VARIANT_ID,Number=.,Type=String,Description="WHO graded variant identifier, one per DRUG, same order">
##INFO=<ID=GENOME_POSITION,Number=1,Type=Integer,Description="Genome position">
"""

        with open(args.header, 'w') as f:
            f.write(header_content)
        print(f"Header file saved to {args.header}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()