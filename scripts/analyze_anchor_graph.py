#!/usr/bin/env python3
"""Analyze and save anchor graph. Compute graph statistics,
characterize nodes, edges, blocks, gaps."""

import os, json, statistics
from collections import defaultdict

def get_genes_ordered(gff_path):
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene': continue
            if 'gene_biotype=protein_coding' not in fields[8]: continue
            name = None
            for part in fields[8].split(';'):
                if part.startswith('gene='):
                    name = part.split('=')[1]
            if not name: continue
            genes.append({
                'name': name,
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6],
                'length': int(fields[4]) - int(fields[3]) + 1,
            })
    genes.sort(key=lambda g: g['start'])
    return genes

def get_all_genes(gff_path):
    """Get ALL protein-coding genes (named + unnamed)."""
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene': continue
            if 'gene_biotype=protein_coding' not in fields[8]: continue
            name = None
            locus = None
            for part in fields[8].split(';'):
                if part.startswith('gene='): name = part.split('=')[1]
                if part.startswith('locus_tag='): locus = part.split('=')[1]
            genes.append({
                'name': name or locus or f"pos_{fields[3]}",
                'has_name': name is not None,
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6],
                'length': int(fields[4]) - int(fields[3]) + 1,
            })
    genes.sort(key=lambda g: g['start'])
    return genes

def main():
    norm_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'normalized')
    out_dir = os.path.join(os.path.dirname(__file__), '..', 'data')

    genomes = [
        "ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
        "pseudomonas_aeruginosa", "saureus_nctc8325",
        "mycobacterium_tuberculosis", "mycoplasma_genitalium",
        "synechocystis_pcc6803", "bacteroides_fragilis",
        "borrelia_burgdorferi"
    ]

    # Load named genes
    named_genes = {}
    all_genes_data = {}
    for name in genomes:
        gff = os.path.join(norm_dir, f"{name}.gff")
        if not os.path.exists(gff): continue
        named_genes[name] = get_genes_ordered(gff)
        all_genes_data[name] = get_all_genes(gff)

    # Build anchor set (≥4 genomes)
    name_counts = defaultdict(int)
    name_genomes = defaultdict(set)
    for gname, genes in named_genes.items():
        for g in genes:
            name_counts[g['name']] += 1
            name_genomes[g['name']].add(gname)

    anchors = {n: list(gs) for n, gs in name_genomes.items() if len(gs) >= 4}

    # Build per-genome anchor chains
    chains = {}
    for gname, genes in named_genes.items():
        chain = [g for g in genes if g['name'] in anchors]
        chains[gname] = chain

    # Build edge set
    edges = defaultdict(lambda: {'genomes': set(), 'gaps': [], 'orientations': []})
    for gname, chain in chains.items():
        for i in range(len(chain) - 1):
            a, b = chain[i], chain[i+1]
            key = (a['name'], b['name'])
            edges[key]['genomes'].add(gname)
            edges[key]['gaps'].append(b['start'] - a['end'])
            edges[key]['orientations'].append((a['strand'], b['strand']))

    # ========== GRAPH STATISTICS ==========
    print(f"{'='*70}")
    print(f"  ANCHOR GRAPH STATISTICS")
    print(f"{'='*70}")

    print(f"\n  NODES (anchor genes):")
    print(f"    Total anchors: {len(anchors)}")
    for min_g in [2, 4, 6, 8, 10]:
        n = sum(1 for gs in name_genomes.values() if len(gs) >= min_g)
        print(f"    In ≥{min_g} genomes: {n}")

    # Node degree distribution
    node_degree = defaultdict(int)
    for (a, b), info in edges.items():
        node_degree[a] += 1
        node_degree[b] += 1
    degrees = list(node_degree.values())
    if degrees:
        print(f"\n    Degree distribution: min={min(degrees)}, median={statistics.median(degrees):.0f}, max={max(degrees)}")
        for d in [1, 2, 3, 4, 5, 10]:
            n = sum(1 for x in degrees if x >= d)
            print(f"      degree ≥{d}: {n} nodes")

    print(f"\n  EDGES (anchor-to-anchor intervals):")
    print(f"    Total unique edges: {len(edges)}")
    conserved = {k: v for k, v in edges.items() if len(v['genomes']) >= 2}
    print(f"    Conserved (≥2 genomes): {len(conserved)}")

    conservation_dist = defaultdict(int)
    for v in edges.values():
        conservation_dist[len(v['genomes'])] += 1
    for n in sorted(conservation_dist.keys()):
        print(f"      in {n} genomes: {conservation_dist[n]}")

    # Gap size distribution for conserved edges
    conserved_gaps = []
    for info in conserved.values():
        conserved_gaps.extend(info['gaps'])
    if conserved_gaps:
        print(f"\n    Gap sizes (conserved edges):")
        print(f"      median={statistics.median(conserved_gaps):.0f}bp")
        print(f"      mean={statistics.mean(conserved_gaps):.0f}bp")
        for bp in [0, 1000, 5000, 10000, 50000]:
            n = sum(1 for g in conserved_gaps if g <= bp)
            print(f"      ≤{bp}bp: {n}/{len(conserved_gaps)} ({n*100//len(conserved_gaps)}%)")

    # Anchor spacing per genome
    print(f"\n  ANCHOR DENSITY per genome:")
    for gname in genomes:
        if gname not in chains: continue
        chain = chains[gname]
        if len(chain) < 2: continue
        spacings = [chain[i+1]['start'] - chain[i]['end'] for i in range(len(chain)-1)]
        genome_len = max(g['end'] for g in all_genes_data.get(gname, [{'end': 1}]))
        total_genes = len(all_genes_data.get(gname, []))
        anchor_frac = len(chain) / max(total_genes, 1)
        print(f"    {gname}: {len(chain)} anchors / {total_genes} genes ({anchor_frac:.0%}), "
              f"median spacing={statistics.median(spacings):.0f}bp")

    # Non-anchor genes between anchors
    print(f"\n  NON-ANCHOR CONTENT (genes between consecutive anchors):")
    for gname in ["ecoli_k12", "salmonella_lt2", "bacillus_subtilis"]:
        if gname not in chains: continue
        chain = chains[gname]
        all_g = all_genes_data.get(gname, [])
        anchor_names = set(g['name'] for g in chain)

        between_counts = []
        for i in range(len(chain) - 1):
            left_end = chain[i]['end']
            right_start = chain[i+1]['start']
            between = [g for g in all_g if g['start'] >= left_end and g['end'] <= right_start
                       and g['name'] not in anchor_names]
            between_counts.append(len(between))

        if between_counts:
            print(f"    {gname}: median {statistics.median(between_counts):.0f} non-anchor genes per interval, "
                  f"max={max(between_counts)}, total intervals={len(between_counts)}")

    # ========== SAVE GRAPH ==========
    graph = {
        'anchors': {name: list(gs) for name, gs in anchors.items()},
        'chains': {},
        'edges': {},
    }

    for gname, chain in chains.items():
        graph['chains'][gname] = [
            {'name': g['name'], 'start': g['start'], 'end': g['end'],
             'strand': g['strand'], 'length': g['length']}
            for g in chain
        ]

    for (a, b), info in edges.items():
        key = f"{a}|{b}"
        graph['edges'][key] = {
            'left': a, 'right': b,
            'genomes': list(info['genomes']),
            'n_genomes': len(info['genomes']),
            'gaps': info['gaps'],
            'orientations': [(s1, s2) for s1, s2 in info['orientations']],
        }

    out_path = os.path.join(out_dir, 'anchor_graph.json')
    with open(out_path, 'w') as f:
        json.dump(graph, f, indent=2)
    print(f"\n  Graph saved: {out_path}")
    print(f"  {len(graph['anchors'])} anchors, {len(graph['edges'])} edges, {len(graph['chains'])} genomes")

if __name__ == "__main__":
    main()
