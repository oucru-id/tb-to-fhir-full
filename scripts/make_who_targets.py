#!/usr/bin/env python3

import argparse
import gzip
import sys
from collections import defaultdict


def load_gene_spans(annotation_table):
    spans = {}

    opener = gzip.open if str(annotation_table).endswith('.gz') else open
    with opener(annotation_table, 'rt') as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 10:
                continue

            chrom, pos, ref, alt, gene, drug = fields[:6]

            if not gene or gene.lower() == 'unknown':
                continue
            if not drug or drug.lower() == 'unknown':
                continue

            try:
                position = int(pos)
            except ValueError:
                continue

            end = position + max(len(ref), 1) - 1

            entry = spans.get(gene)
            if entry is None:
                spans[gene] = {
                    'chrom': chrom,
                    'start': position,
                    'end': end,
                    'drugs': set(),
                }
                entry = spans[gene]
            else:
                entry['start'] = min(entry['start'], position)
                entry['end'] = max(entry['end'], end)

            for one_drug in drug.split(','):
                one_drug = one_drug.strip()
                if one_drug and one_drug.lower() != 'unknown':
                    entry['drugs'].add(one_drug)

    return spans


def load_mask(path):
    if not path:
        return []

    intervals = []
    try:
        with open(path) as handle:
            for line in handle:
                if line.startswith(('#', 'track', 'browser')):
                    continue
                fields = line.split()
                if len(fields) < 3:
                    continue
                try:
                    intervals.append((fields[0], int(fields[1]), int(fields[2])))
                except ValueError:
                    continue
    except FileNotFoundError:
        print(f"Warning: repetitive regions file not found: {path}", file=sys.stderr)
        return []

    return intervals


def subtract(chrom, start, end, mask):
    pieces = [(start, end)]

    for m_chrom, m_start, m_end in mask:
        if m_chrom != chrom:
            continue
        survivors = []
        for p_start, p_end in pieces:
            if m_end <= p_start or m_start >= p_end:
                survivors.append((p_start, p_end))
                continue
            if m_start > p_start:
                survivors.append((p_start, m_start))
            if m_end < p_end:
                survivors.append((m_end, p_end))
        pieces = survivors
        if not pieces:
            break

    return pieces


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--annotation-table', required=True,
                        help='WHO catalogue annotation table (.tsv or .tsv.gz)')
    parser.add_argument('--output', required=True, help='Output BED path')
    parser.add_argument('--repetitive-regions', default=None,
                        help='BED of regions to exclude (pe/ppe, IS elements)')
    parser.add_argument('--promoter-padding', type=int, default=200,
                        help='Bases to extend upstream/downstream of each gene span '
                             '(default: 200, to capture promoter variants)')
    args = parser.parse_args()

    spans = load_gene_spans(args.annotation_table)
    if not spans:
        print('ERROR: no gene spans found in the annotation table', file=sys.stderr)
        sys.exit(1)

    mask = load_mask(args.repetitive_regions)

    rows = []
    dropped_entirely = []

    for gene, entry in spans.items():
        chrom = entry['chrom']
        start = max(0, entry['start'] - 1 - args.promoter_padding)
        end = entry['end'] + args.promoter_padding

        drugs = ';'.join(sorted(entry['drugs'])) or 'unknown'
        pieces = subtract(chrom, start, end, mask)

        if not pieces:
            dropped_entirely.append(gene)
            continue

        for p_start, p_end in pieces:
            if p_end > p_start:
                rows.append((chrom, p_start, p_end, gene, drugs))

    rows.sort(key=lambda r: (r[0], r[1], r[2]))

    with open(args.output, 'w') as out:
        out.write('#chrom\tstart\tend\tgene\tdrugs\n')
        for chrom, start, end, gene, drugs in rows:
            out.write(f'{chrom}\t{start}\t{end}\t{gene}\t{drugs}\n')

    drug_genes = defaultdict(set)
    for _, _, _, gene, drugs in rows:
        for drug in drugs.split(';'):
            drug_genes[drug].add(gene)

    print(f'Wrote {len(rows)} intervals covering {len(spans) - len(dropped_entirely)} genes '
          f'to {args.output}')
    print(f'Drugs represented: {len(drug_genes)}')
    for drug in sorted(drug_genes):
        genes = sorted(drug_genes[drug])
        preview = ', '.join(genes[:6]) + (' ...' if len(genes) > 6 else '')
        print(f'  {drug:<16} {len(genes):>3} loci  ({preview})')

    if dropped_entirely:
        print(f'\nWarning: {len(dropped_entirely)} genes fell entirely inside the '
              f'repetitive-region mask and have no target interval:', file=sys.stderr)
        print('  ' + ', '.join(sorted(dropped_entirely)), file=sys.stderr)
        print('  Drugs relying solely on these loci can never be assessable.',
              file=sys.stderr)


if __name__ == '__main__':
    main()
