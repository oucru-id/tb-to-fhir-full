#!/usr/bin/env python3

import argparse
import gzip
import json
import sys
from collections import defaultdict


def _open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def load_targets(path):
    intervals = {}
    gene_drugs = defaultdict(set)
    gene_gating = defaultdict(set)
    has_gating = False

    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 5:
                continue
            chrom, start, end, gene, drugs = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4]
            drug_list = [d.strip() for d in drugs.split(';') if d.strip() and d.strip().lower() != 'unknown']
            intervals[(chrom, start, end)] = (gene, drug_list)
            gene_drugs[gene].update(drug_list)
            if len(fields) > 5:
                has_gating = True
                gene_gating[gene].update(d.strip() for d in fields[5].split(';') if d.strip() not in ('', '.'))

    return intervals, gene_drugs, (gene_gating if has_gating else None)


def load_regions(path):
    depths = {}
    with _open(path) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4:
                continue
            try:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                mean_depth = float(fields[-1])
            except ValueError:
                continue
            depths[(chrom, start, end)] = mean_depth
    return depths


def load_positions(path):
    positions = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 5:
                continue
            key = (fields[0], int(fields[1]), int(fields[2]))
            positions[key] = (fields[3], [d for d in fields[4].split(';') if d])
    return positions


def describe(start, end):
    return str(end) if end - start == 1 else f'{start + 1}-{end}'


def load_thresholds(path, min_depth):

    breadths = {}
    target_col = None

    with _open(path) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if not fields:
                continue

            if line.startswith('#'):
                wanted = f'{min_depth}X'
                for idx, name in enumerate(fields):
                    if name.strip().upper() == wanted.upper():
                        target_col = idx
                        break
                continue

            if target_col is None or len(fields) <= target_col:
                continue

            try:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                bases_at_depth = int(fields[target_col])
            except ValueError:
                continue

            span = end - start
            breadths[(chrom, start, end)] = (bases_at_depth / span) if span > 0 else 0.0

    return breadths


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sample-id', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--targets')
    parser.add_argument('--regions')
    parser.add_argument('--thresholds')
    parser.add_argument('--min-depth', type=int, default=10)
    parser.add_argument('--min-breadth', type=float, default=0.95)
    parser.add_argument('--positions', default=None,
                        help='BED of WHO grade 1-2 positions (make_who_targets.py --positions-output)')
    parser.add_argument('--position-thresholds', default=None,
                        help='mosdepth thresholds.bed.gz computed over --positions')
    parser.add_argument('--not-assessed', default=None,
                        help='Emit an explicit "coverage not assessed" summary with this reason')
    args = parser.parse_args()

    if args.not_assessed or not (args.targets and args.regions and args.thresholds):
        reason = args.not_assessed or 'coverage inputs not provided'
        summary = {
            'sample_id': args.sample_id,
            'coverage_assessed': False,
            'reason': reason,
            'min_depth': args.min_depth,
            'min_breadth': args.min_breadth,
            'genes': {},
            'drugs': {},
        }
        with open(args.output, 'w') as out:
            json.dump(summary, out, indent=2)
        print(f'Coverage not assessed for {args.sample_id}: {reason}')
        return

    intervals, gene_drugs, gene_gating = load_targets(args.targets)
    depths = load_regions(args.regions)
    breadths = load_thresholds(args.thresholds, args.min_depth)

    if not breadths:
        print(f'ERROR: no usable threshold columns for {args.min_depth}X in '
              f'{args.thresholds}', file=sys.stderr)
        sys.exit(1)

    gene_bases = defaultdict(int)
    gene_covered = defaultdict(float)
    gene_depth_weighted = defaultdict(float)

    for key, (gene, _drugs) in intervals.items():
        chrom, start, end = key
        span = end - start
        if span <= 0:
            continue
        breadth = breadths.get(key)
        if breadth is None:
            breadth = 0.0
        gene_bases[gene] += span
        gene_covered[gene] += breadth * span
        gene_depth_weighted[gene] += depths.get(key, 0.0) * span

    genes = {}
    for gene, bases in gene_bases.items():
        breadth = (gene_covered[gene] / bases) if bases else 0.0
        mean_depth = (gene_depth_weighted[gene] / bases) if bases else 0.0
        genes[gene] = {
            'bases': bases,
            'mean_depth': round(mean_depth, 2),
            f'breadth_at_{args.min_depth}x': round(breadth, 4),
            'adequate': breadth >= args.min_breadth,
            'drugs': sorted(gene_drugs.get(gene, [])),
        }

    use_positions = bool(args.positions and args.position_thresholds and gene_gating is not None)
    pair_positions = defaultdict(list)
    if use_positions:
        position_breadths = load_thresholds(args.position_thresholds, args.min_depth)
        for key, (gene, position_drugs) in load_positions(args.positions).items():
            covered = position_breadths.get(key, 0.0) >= 1.0
            for drug in position_drugs:
                pair_positions[(drug, gene)].append((key[1], key[2], covered))
            if gene in genes:
                genes[gene].setdefault('graded_positions', 0)
                genes[gene].setdefault('graded_positions_at_min_depth', 0)
                genes[gene]['graded_positions'] += 1
                genes[gene]['graded_positions_at_min_depth'] += int(covered)

    drug_to_genes = defaultdict(set)
    for gene, drugs in gene_drugs.items():
        for drug in drugs:
            drug_to_genes[drug].add(gene)

    drugs = {}
    for drug, drug_genes in drug_to_genes.items():
        if gene_gating is None:
            gating = set(drug_genes)
        else:
            gating = {g for g in drug_genes if drug in gene_gating.get(g, set())}
        positional = {g for g in gating if (drug, g) in pair_positions}
        uncovered = {g: [describe(s, e) for s, e, ok in sorted(pair_positions[(drug, g)]) if not ok]
                     for g in sorted(positional)}
        inadequate = sorted([g for g in gating - positional if not genes.get(g, {}).get('adequate', False)]
                            + [g for g in positional if uncovered[g]])
        missing = sorted(g for g in gating - positional if g not in genes)
        other = sorted(g for g in drug_genes - gating if not genes.get(g, {}).get('adequate', False))
        drugs[drug] = {
            'loci': sorted(drug_genes),
            'gating_loci': sorted(gating),
            'assessable': bool(gating) and not inadequate and not missing,
            'inadequate_loci': inadequate,
            'loci_without_coverage_data': missing,
            'inadequate_other_loci': other,
        }
        if use_positions:
            drugs[drug]['position_gated_loci'] = sorted(positional)
            drugs[drug]['uncovered_positions'] = {g: p for g, p in uncovered.items() if p}
        if not gating:
            drugs[drug]['reason'] = ('no locus carries a WHO catalogue variant graded as associated '
                                     'with resistance to this drug')

    if use_positions:
        gating_rule = ('WHO grade 1-2 positions for loci without graded loss-of-function variants; '
                       'whole locus otherwise')
    elif gene_gating is not None:
        gating_rule = 'loci with WHO grade 1-2 variants'
    else:
        gating_rule = 'all catalogue loci'
    summary = {
        'sample_id': args.sample_id,
        'coverage_assessed': True,
        'min_depth': args.min_depth,
        'min_breadth': args.min_breadth,
        'gating': gating_rule,
        'genes': genes,
        'drugs': drugs,
    }

    with open(args.output, 'w') as out:
        json.dump(summary, out, indent=2)

    assessable = sum(1 for d in drugs.values() if d['assessable'])
    print(f'{args.sample_id}: {len(genes)} loci summarized; '
          f'{assessable}/{len(drugs)} drugs assessable at '
          f'{args.min_breadth:.0%} breadth @ {args.min_depth}x')


if __name__ == '__main__':
    main()
