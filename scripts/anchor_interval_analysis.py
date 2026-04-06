#!/usr/bin/env python3
"""Analyze anchor-bounded intervals: expected vs observed gene counts.
Find intervals where we predict too few or too many genes."""

import os, json, statistics
from collections import defaultdict

def get_all_genes(gff_path):
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene': continue
            if 'gene_biotype=protein_coding' not in fields[8]: continue
            name = None
            for part in fields[8].split(';'):
                if part.startswith('gene='): name = part.split('=')[1]
            genes.append({
                'name': name,
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6],
            })
    genes.sort(key=lambda g: g['start'])
    return genes

def count_genes_between(all_genes, left_end, right_start):
    """Count genes strictly between two positions."""
    return sum(1 for g in all_genes if g['start'] > left_end and g['end'] < right_start)

def main():
    norm_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'normalized')
    graph_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'anchor_graph.json')

    with open(graph_path) as f:
        graph = json.load(f)

    genomes = [
        "ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
        "pseudomonas_aeruginosa", "mycobacterium_tuberculosis",
        "synechocystis_pcc6803", "bacteroides_fragilis",
    ]

    # Load reference genes per genome
    ref_genes = {}
    for gname in genomes:
        gff = os.path.join(norm_dir, f"{gname}.gff")
        if os.path.exists(gff):
            ref_genes[gname] = get_all_genes(gff)

    # For each conserved edge (anchor pair in ≥3 genomes),
    # compute gene count in each genome's interval
    edges = graph['edges']
    chains = graph['chains']

    # Build per-genome anchor positions
    anchor_pos = {}  # genome -> {anchor_name: (start, end)}
    for gname, chain in chains.items():
        anchor_pos[gname] = {}
        for g in chain:
            anchor_pos[gname][g['name']] = (g['start'], g['end'])

    interval_data = []  # (left, right, genome, n_genes, span_bp)

    for edge_key, info in edges.items():
        if info['n_genomes'] < 3: continue
        left, right = info['left'], info['right']

        gene_counts = []
        spans = []

        for gname in info['genomes']:
            if gname not in ref_genes: continue
            if gname not in anchor_pos: continue
            if left not in anchor_pos[gname] or right not in anchor_pos[gname]:
                continue

            l_start, l_end = anchor_pos[gname][left]
            r_start, r_end = anchor_pos[gname][right]

            if l_end >= r_start: continue  # overlapping anchors

            n_between = count_genes_between(ref_genes[gname], l_end, r_start)
            span = r_start - l_end

            gene_counts.append(n_between)
            spans.append(span)

            interval_data.append({
                'left': left, 'right': right,
                'genome': gname, 'n_genes': n_between, 'span': span,
            })

    print(f"Analyzed {len(interval_data)} anchor-bounded intervals")
    print(f"From {sum(1 for e in edges.values() if e['n_genomes'] >= 3)} conserved edges (≥3 genomes)")

    # Group by edge: compare gene counts across genomes
    edge_groups = defaultdict(list)
    for iv in interval_data:
        edge_groups[(iv['left'], iv['right'])].append(iv)

    # Find consistent intervals (low variance) vs variable ones
    consistent = 0
    variable = 0
    zero_gene = 0
    multi_gene = 0

    expected_gene_counts = {}  # (left, right) -> median gene count

    for (left, right), intervals in edge_groups.items():
        counts = [iv['n_genes'] for iv in intervals]
        if len(counts) < 2: continue

        med = statistics.median(counts)
        expected_gene_counts[(left, right)] = med

        if max(counts) - min(counts) <= 1:
            consistent += 1
        else:
            variable += 1

        if med == 0:
            zero_gene += 1
        elif med >= 3:
            multi_gene += 1

    print(f"\nInterval consistency:")
    print(f"  Consistent (range ≤1): {consistent}")
    print(f"  Variable (range >1): {variable}")
    print(f"\nGene count distribution:")
    print(f"  0 genes between anchors: {zero_gene} intervals (adjacent operonic genes)")
    print(f"  ≥3 genes between anchors: {multi_gene} intervals")

    # Distribution of expected gene counts
    all_medians = list(expected_gene_counts.values())
    if all_medians:
        print(f"\nExpected gene count per interval:")
        for n in [0, 1, 2, 3, 5, 10, 20, 50]:
            k = sum(1 for m in all_medians if m <= n)
            print(f"  ≤{n} genes: {k}/{len(all_medians)} ({k*100//len(all_medians)}%)")

    # KEY OUTPUT: save expected gene counts for use in annotation
    out_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'anchor_intervals.json')
    intervals_json = {}
    for (left, right), med in expected_gene_counts.items():
        key = f"{left}|{right}"
        group = edge_groups[(left, right)]
        intervals_json[key] = {
            'left': left,
            'right': right,
            'expected_genes': med,
            'n_genomes': len(group),
            'gene_counts': [iv['n_genes'] for iv in group],
            'spans': [iv['span'] for iv in group],
        }

    with open(out_path, 'w') as f:
        json.dump(intervals_json, f, indent=2)
    print(f"\nSaved {len(intervals_json)} interval expectations to {out_path}")

if __name__ == "__main__":
    main()
