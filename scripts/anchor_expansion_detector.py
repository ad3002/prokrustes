#!/usr/bin/env python3
"""Detect insertions between conserved anchors.
For anchor pairs that are normally adjacent (0-2 genes),
find genomes where the interval is strongly expanded."""

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

def find_known_phage_islands(gff_path):
    phage_genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS': continue
            attrs = fields[8].lower()
            if any(kw in attrs for kw in ['phage', 'prophage', 'capsid', 'terminase',
                    'integrase', 'portal', 'tail', 'baseplate']):
                phage_genes.append((int(fields[3]), int(fields[4])))
    phage_genes.sort()
    clusters = []
    cur = [phage_genes[0]] if phage_genes else []
    for g in phage_genes[1:]:
        if g[0] - cur[-1][1] < 15000:
            cur.append(g)
        else:
            if len(cur) >= 2:
                clusters.append((cur[0][0], cur[-1][1]))
            cur = [g]
    if len(cur) >= 2:
        clusters.append((cur[0][0], cur[-1][1]))
    return clusters

def in_island(s, e, islands):
    for ps, pe in islands:
        if s < pe and e > ps: return True
    return False

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"
    norm_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'normalized')
    graph_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'anchor_graph.json')

    with open(graph_path) as f:
        graph = json.load(f)

    genomes = [
        "ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
        "pseudomonas_aeruginosa", "saureus_nctc8325",
        "mycobacterium_tuberculosis", "synechocystis_pcc6803",
        "bacteroides_fragilis", "borrelia_burgdorferi",
    ]

    # Load normalized genes
    norm_genes = {}
    for gname in genomes:
        gff = os.path.join(norm_dir, f"{gname}.gff")
        if os.path.exists(gff):
            norm_genes[gname] = get_all_genes(gff)

    # Load known phage islands (from original GFF)
    known_phage = {}
    for gname in genomes:
        gff = os.path.join(data_dir, f"{gname}.gff")
        if os.path.exists(gff):
            known_phage[gname] = find_known_phage_islands(gff)

    # Build anchor positions per genome
    anchor_pos = {}
    for gname, chain in graph['chains'].items():
        anchor_pos[gname] = {g['name']: (g['start'], g['end']) for g in chain}

    # For each conserved edge, compute gene count per genome
    edges = graph['edges']
    expansions = []

    for edge_key, info in edges.items():
        if info['n_genomes'] < 3: continue
        left, right = info['left'], info['right']

        # Collect gene counts per genome
        genome_counts = {}
        for gname in info['genomes']:
            if gname not in norm_genes or gname not in anchor_pos: continue
            if left not in anchor_pos[gname] or right not in anchor_pos[gname]: continue

            l_start, l_end = anchor_pos[gname][left]
            r_start, r_end = anchor_pos[gname][right]
            if l_end >= r_start: continue

            n_between = sum(1 for g in norm_genes[gname]
                          if g['start'] > l_end and g['end'] < r_start)
            span = r_start - l_end
            genome_counts[gname] = {'n_genes': n_between, 'span': span}

        if len(genome_counts) < 3: continue

        counts = [v['n_genes'] for v in genome_counts.values()]
        median_count = statistics.median(counts)
        spans = [v['span'] for v in genome_counts.values()]
        median_span = statistics.median(spans)

        # Find expanded genomes: significantly more genes than expected
        for gname, data in genome_counts.items():
            excess_genes = data['n_genes'] - median_count
            excess_span = data['span'] - median_span

            # Expansion: ≥5 more genes than median, or ≥5000bp more span
            if excess_genes >= 5 or (excess_span >= 5000 and data['n_genes'] >= 3):
                # Check if this overlaps known phage
                l_end = anchor_pos[gname][left][1]
                r_start = anchor_pos[gname][right][0]
                is_phage = in_island(l_end, r_start, known_phage.get(gname, []))

                expansions.append({
                    'left': left, 'right': right,
                    'genome': gname,
                    'expected_genes': median_count,
                    'observed_genes': data['n_genes'],
                    'excess_genes': excess_genes,
                    'expected_span': median_span,
                    'observed_span': data['span'],
                    'excess_span': excess_span,
                    'is_known_phage': is_phage,
                    'left_pos': anchor_pos[gname][left][1],
                    'right_pos': anchor_pos[gname][right][0],
                })

    print(f"{'='*70}")
    print(f"  ANCHOR EXPANSION DETECTOR")
    print(f"{'='*70}")
    print(f"\n  Total expansions detected: {len(expansions)}")

    # Group by genome
    by_genome = defaultdict(list)
    for e in expansions:
        by_genome[e['genome']].append(e)

    for gname in genomes:
        exps = by_genome.get(gname, [])
        if not exps: continue
        n_phage = sum(1 for e in exps if e['is_known_phage'])
        n_other = len(exps) - n_phage
        total_excess_bp = sum(e['excess_span'] for e in exps)
        print(f"\n  {gname}: {len(exps)} expansions ({n_phage} phage, {n_other} other), {total_excess_bp//1000}kb excess")
        for e in sorted(exps, key=lambda x: -x['excess_genes'])[:5]:
            label = "PHAGE" if e['is_known_phage'] else "other"
            print(f"    {e['left']:>10}..{e['right']:<10} +{e['excess_genes']} genes "
                  f"(expected {e['expected_genes']:.0f}, got {e['observed_genes']}) "
                  f"+{e['excess_span']//1000}kb  [{label}]")

    # Summary: how well do expansions correlate with known phage?
    total_phage = sum(1 for e in expansions if e['is_known_phage'])
    total_other = len(expansions) - total_phage
    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"  Total expansions: {len(expansions)}")
    print(f"  Known phage:      {total_phage} ({total_phage*100//max(len(expansions),1)}%)")
    print(f"  Other (IS/acc.):  {total_other} ({total_other*100//max(len(expansions),1)}%)")

    # Save expansions
    out_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'anchor_expansions.json')
    with open(out_path, 'w') as f:
        json.dump(expansions, f, indent=2)
    print(f"\n  Saved: {out_path}")

if __name__ == "__main__":
    main()
