#!/usr/bin/env python3
"""Measure local codon usage variation across anchor-defined blocks.
Key question: is there enough variation WITHIN a genome to justify
local coding models instead of one global model?"""

import os, json, statistics, math
from collections import defaultdict

def get_all_genes_with_seq(gff_path, fasta_path):
    """Get protein-coding genes with DNA sequences."""
    seqs = {}
    name = ""
    parts = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if name: seqs[name] = ''.join(parts)
                name = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip().upper())
    if name: seqs[name] = ''.join(parts)

    genome = list(seqs.values())[0]
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

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

            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            dna = genome[start:end]
            if strand == '-':
                dna = ''.join(comp.get(c, 'N') for c in reversed(dna))

            genes.append({
                'name': name,
                'start': start + 1, 'end': end, 'strand': strand,
                'dna': dna, 'length': len(dna),
            })
    genes.sort(key=lambda g: g['start'])
    return genes

def hex_profile(dna_seqs):
    """Build hexamer frequency profile from DNA sequences."""
    counts = defaultdict(int)
    total = 0
    nt = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for seq in dna_seqs:
        for i in range(0, len(seq) - 5, 3):  # in-frame hexamers
            h = seq[i:i+6]
            if all(c in nt for c in h):
                counts[h] += 1
                total += 1
    # Normalize
    profile = {}
    for h, c in counts.items():
        profile[h] = c / max(total, 1)
    return profile

def kl_divergence(p, q):
    """KL(P || Q) for two profiles."""
    kl = 0.0
    all_keys = set(list(p.keys()) + list(q.keys()))
    n = len(all_keys)
    for k in all_keys:
        pk = p.get(k, 0.5/n)
        qk = q.get(k, 0.5/n)
        if pk > 0 and qk > 0:
            kl += pk * math.log(pk / qk)
    return kl

def main():
    norm_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'normalized')
    graph_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'anchor_graph.json')

    with open(graph_path) as f:
        graph = json.load(f)

    genomes_to_test = ["ecoli_k12", "pseudomonas_aeruginosa", "mycobacterium_tuberculosis",
                       "bacillus_subtilis", "salmonella_lt2"]

    for gname in genomes_to_test:
        gff = os.path.join(norm_dir, f"{gname}.gff")
        fasta = os.path.join(norm_dir, f"{gname}.fasta")
        if not os.path.exists(gff): continue

        genes = get_all_genes_with_seq(gff, fasta)
        if not genes: continue

        # Build anchor positions
        chain = graph['chains'].get(gname, [])
        anchor_names = set(g['name'] for g in chain)
        anchor_positions = [(g['start'], g['end'], g['name']) for g in chain]
        anchor_positions.sort()

        # Split genome into blocks between SPARSE anchors (skip adjacent ones)
        # Take every 10th anchor to get ~10kb blocks
        sparse_anchors = anchor_positions[::10]
        if len(sparse_anchors) < 5:
            sparse_anchors = anchor_positions[::5]

        blocks = []
        for i in range(len(sparse_anchors) - 1):
            left_end = sparse_anchors[i][1]
            right_start = sparse_anchors[i+1][0]
            block_genes = [g for g in genes if g['start'] >= left_end and g['end'] <= right_start
                          and g['length'] >= 300]
            if len(block_genes) >= 5:
                blocks.append({
                    'left': sparse_anchors[i][2],
                    'right': sparse_anchors[i+1][2],
                    'genes': block_genes,
                    'n_genes': len(block_genes),
                    'start': left_end,
                    'end': right_start,
                })

        if len(blocks) < 3: continue

        # Build global hex profile
        all_dna = [g['dna'] for g in genes if g['length'] >= 300]
        global_profile = hex_profile(all_dna)

        # Build per-block hex profiles and compare to global
        block_kls = []
        for b in blocks:
            block_dna = [g['dna'] for g in b['genes']]
            block_prof = hex_profile(block_dna)
            kl = kl_divergence(block_prof, global_profile)
            block_kls.append(kl)
            b['kl_from_global'] = kl

        # Also compute GC per block
        for b in blocks:
            all_seq = ''.join(g['dna'] for g in b['genes'])
            gc = sum(1 for c in all_seq if c in 'GC') / max(len(all_seq), 1)
            b['gc'] = gc

        genome_gc = sum(1 for c in ''.join(all_dna) if c in 'GC') / max(sum(len(s) for s in all_dna), 1)

        print(f"\n{'='*60}")
        print(f"  {gname}: {len(blocks)} blocks, {len(genes)} genes")
        print(f"{'='*60}")
        print(f"  Genome GC: {genome_gc:.3f}")
        print(f"  Block KL from global hex: median={statistics.median(block_kls):.4f}, "
              f"max={max(block_kls):.4f}, min={min(block_kls):.4f}")

        # Identify anomalous blocks (high KL = different codon usage)
        kl_sorted = sorted(block_kls)
        kl_p90 = kl_sorted[len(kl_sorted) * 90 // 100]
        anomalous = [b for b in blocks if b['kl_from_global'] > kl_p90]

        print(f"  Anomalous blocks (top 10% KL): {len(anomalous)}")
        for b in sorted(anomalous, key=lambda x: -x['kl_from_global'])[:5]:
            print(f"    {b['left']:>10}..{b['right']:<10}: {b['n_genes']} genes, "
                  f"KL={b['kl_from_global']:.4f}, GC={b['gc']:.3f} (genome={genome_gc:.3f})")

        # Key question: do anomalous blocks correspond to phage/IS regions?
        # Check GC shift
        gc_shifts = [abs(b['gc'] - genome_gc) for b in blocks]
        anomalous_gc_shifts = [abs(b['gc'] - genome_gc) for b in anomalous]
        normal_gc_shifts = [abs(b['gc'] - genome_gc) for b in blocks if b not in anomalous]

        if anomalous_gc_shifts and normal_gc_shifts:
            print(f"\n  GC shift from genome mean:")
            print(f"    Normal blocks: median |ΔGC| = {statistics.median(normal_gc_shifts):.4f}")
            print(f"    Anomalous blocks: median |ΔGC| = {statistics.median(anomalous_gc_shifts):.4f}")

        # Would local models help?
        # If anomalous blocks have KL > 0.01, that's ~5% codon usage difference
        # which could affect hex scoring significantly
        significant = sum(1 for kl in block_kls if kl > 0.01)
        print(f"\n  Blocks with KL > 0.01 (significant codon shift): {significant}/{len(blocks)}")

if __name__ == "__main__":
    main()
