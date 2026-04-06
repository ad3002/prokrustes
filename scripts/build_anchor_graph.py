#!/usr/bin/env python3
"""Build anchor-based genome segmentation and synteny graph.
Anchors = shared gene names (orthologs by name).
Edges = observed anchor-to-anchor intervals in each genome."""

import os
from collections import defaultdict

def get_genes_ordered(gff_path):
    """Get ordered list of protein-coding genes with names, positions, strands."""
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
            })
    genes.sort(key=lambda g: g['start'])
    return genes

def build_anchor_chain(genes, anchor_set):
    """Extract ordered chain of anchor genes from a genome."""
    chain = []
    for g in genes:
        if g['name'] in anchor_set:
            chain.append(g)
    return chain

def build_anchor_pairs(chain):
    """Build list of consecutive anchor pairs with interval info."""
    pairs = []
    for i in range(len(chain) - 1):
        a = chain[i]
        b = chain[i + 1]
        gap = b['start'] - a['end']
        pairs.append({
            'left': a['name'],
            'right': b['name'],
            'left_strand': a['strand'],
            'right_strand': b['strand'],
            'gap_bp': gap,
            'left_pos': a['start'],
            'right_pos': b['start'],
        })
    return pairs

def main():
    norm_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'normalized')

    genomes = [
        ("ecoli_k12", "E. coli"),
        ("salmonella_lt2", "Salmonella"),
        ("bacillus_subtilis", "B. subtilis"),
        ("pseudomonas_aeruginosa", "Pseudomonas"),
        ("saureus_nctc8325", "S. aureus"),
        ("mycobacterium_tuberculosis", "M. tuberculosis"),
        ("mycoplasma_genitalium", "Mycoplasma"),
        ("synechocystis_pcc6803", "Synechocystis"),
        ("bacteroides_fragilis", "Bacteroides"),
        ("borrelia_burgdorferi", "Borrelia"),
    ]

    # Load all genomes
    all_genes = {}
    for name, label in genomes:
        gff = os.path.join(norm_dir, f"{name}.gff")
        if not os.path.exists(gff): continue
        genes = get_genes_ordered(gff)
        all_genes[name] = genes
        print(f"{label}: {len(genes)} named genes")

    # Find shared gene names across all genomes
    name_sets = {}
    for name in all_genes:
        name_sets[name] = set(g['name'] for g in all_genes[name])

    # Core anchors: present in ≥6 genomes
    name_counts = defaultdict(int)
    for ns in name_sets.values():
        for n in ns:
            name_counts[n] += 1

    for min_genomes in [2, 4, 6, 8, 10]:
        n = sum(1 for c in name_counts.values() if c >= min_genomes)
        print(f"\nAnchors in ≥{min_genomes} genomes: {n}")

    # Use anchors present in ≥4 genomes
    anchors = set(n for n, c in name_counts.items() if c >= 4)
    print(f"\nUsing {len(anchors)} anchors (≥4 genomes)")

    # Build anchor chains per genome
    chains = {}
    for name in all_genes:
        chain = build_anchor_chain(all_genes[name], anchors)
        chains[name] = chain
        print(f"  {name}: {len(chain)} anchors in chain")

    # Build anchor pair graph
    # Edge = (left_anchor, right_anchor) → set of genomes where this adjacency exists
    edge_genomes = defaultdict(set)
    edge_info = defaultdict(list)

    for name, chain in chains.items():
        pairs = build_anchor_pairs(chain)
        for p in pairs:
            edge_key = (p['left'], p['right'])
            edge_genomes[edge_key].add(name)
            edge_info[edge_key].append({
                'genome': name,
                'gap_bp': p['gap_bp'],
                'left_strand': p['left_strand'],
                'right_strand': p['right_strand'],
            })

    # Statistics
    total_edges = len(edge_genomes)
    conserved_edges = {k: v for k, v in edge_genomes.items() if len(v) >= 2}
    highly_conserved = {k: v for k, v in edge_genomes.items() if len(v) >= 6}

    print(f"\n{'='*60}")
    print(f"  ANCHOR GRAPH STATISTICS")
    print(f"{'='*60}")
    print(f"  Total anchor pairs (edges): {total_edges}")
    print(f"  Conserved (≥2 genomes):     {len(conserved_edges)}")
    print(f"  Highly conserved (≥6):      {len(highly_conserved)}")

    # Edge conservation distribution
    conservation_dist = defaultdict(int)
    for v in edge_genomes.values():
        conservation_dist[len(v)] += 1
    print(f"\n  Conservation distribution:")
    for n in sorted(conservation_dist.keys()):
        print(f"    in {n} genomes: {conservation_dist[n]} edges")

    # Show top conserved edges
    print(f"\n  Top conserved anchor pairs (≥8 genomes):")
    for (left, right), gset in sorted(edge_genomes.items(), key=lambda x: -len(x[1])):
        if len(gset) < 8: continue
        infos = edge_info[(left, right)]
        gaps = [i['gap_bp'] for i in infos]
        import statistics
        med_gap = statistics.median(gaps)
        print(f"    {left:>10} → {right:<10}: {len(gset)} genomes, median gap={med_gap:.0f}bp")

    # Synteny blocks: find runs of conserved anchor pairs (same order in multiple genomes)
    print(f"\n{'='*60}")
    print(f"  SYNTENY BLOCKS (E. coli vs Salmonella)")
    print(f"{'='*60}")

    ec_chain = chains.get('ecoli_k12', [])
    sal_chain = chains.get('salmonella_lt2', [])

    if ec_chain and sal_chain:
        ec_names = [g['name'] for g in ec_chain]
        sal_names = [g['name'] for g in sal_chain]
        sal_idx = {name: i for i, name in enumerate(sal_names)}

        # Find collinear runs
        blocks = []
        current_block = []
        last_sal_i = -1

        for ec_i, ec_name in enumerate(ec_names):
            if ec_name in sal_idx:
                si = sal_idx[ec_name]
                if last_sal_i >= 0 and si == last_sal_i + 1:
                    # Collinear continuation
                    current_block.append(ec_name)
                else:
                    # Break
                    if len(current_block) >= 3:
                        blocks.append(current_block)
                    current_block = [ec_name]
                last_sal_i = si
            else:
                if len(current_block) >= 3:
                    blocks.append(current_block)
                current_block = []
                last_sal_i = -1

        if len(current_block) >= 3:
            blocks.append(current_block)

        blocks.sort(key=lambda b: -len(b))
        total_anchors_in_blocks = sum(len(b) for b in blocks)

        print(f"  Collinear blocks (≥3 anchors): {len(blocks)}")
        print(f"  Total anchors in blocks: {total_anchors_in_blocks}/{len(ec_chain)}")
        print(f"\n  Largest blocks:")
        for b in blocks[:15]:
            print(f"    {len(b)} anchors: {b[0]} → ... → {b[-1]}")

        # Inversions: runs in reverse order
        inv_blocks = []
        current_inv = []
        last_sal_i = -1

        for ec_name in ec_names:
            if ec_name in sal_idx:
                si = sal_idx[ec_name]
                if last_sal_i >= 0 and si == last_sal_i - 1:
                    current_inv.append(ec_name)
                else:
                    if len(current_inv) >= 3:
                        inv_blocks.append(current_inv)
                    current_inv = [ec_name]
                last_sal_i = si

        if len(current_inv) >= 3:
            inv_blocks.append(current_inv)

        print(f"\n  Inverted blocks (≥3 anchors): {len(inv_blocks)}")
        for b in inv_blocks[:5]:
            print(f"    {len(b)} anchors: {b[0]} → ... → {b[-1]}")

if __name__ == "__main__":
    main()
