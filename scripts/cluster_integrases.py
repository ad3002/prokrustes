#!/usr/bin/env python3
"""Cluster integrases by protein similarity, build family-specific models,
test separation quality per family."""

import os, math, statistics
from collections import defaultdict

def read_fasta(path):
    seqs = {}
    name = ""
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name: seqs[name] = ''.join(parts)
                name = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip().upper())
    if name: seqs[name] = ''.join(parts)
    return seqs

aa_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def translate(dna):
    return ''.join(aa_table.get(dna[i:i+3], 'X') for i in range(0, len(dna) - 2, 3))

def protein_3mer_profile(proteins):
    """Build amino acid 3-mer log-odds vs uniform."""
    counts = defaultdict(int)
    total = 0
    for prot in proteins:
        for i in range(len(prot) - 2):
            tri = prot[i:i+3]
            if '*' in tri or 'X' in tri: continue
            counts[tri] += 1
            total += 1
    n_possible = 20**3  # 8000
    bg = 1.0 / n_possible
    profile = {}
    for tri, cnt in counts.items():
        freq = (cnt + 0.1) / (total + n_possible * 0.1)
        profile[tri] = math.log(freq / bg)
    return profile

def score_protein(prot, profile):
    total = 0.0
    n = 0
    for i in range(len(prot) - 2):
        tri = prot[i:i+3]
        if tri in profile:
            total += profile[tri]
            n += 1
    return total / n if n > 0 else 0.0

def protein_similarity(p1, p2):
    """Simple protein similarity: shared 3-mer fraction."""
    s1 = set(p1[i:i+3] for i in range(len(p1) - 2))
    s2 = set(p2[i:i+3] for i in range(len(p2) - 2))
    if not s1 or not s2: return 0
    return len(s1 & s2) / min(len(s1), len(s2))

def cluster_proteins(proteins, names, threshold=0.35):
    """Simple single-linkage clustering by 3-mer similarity."""
    n = len(proteins)
    clusters = list(range(n))  # each in own cluster

    for i in range(n):
        for j in range(i + 1, n):
            sim = protein_similarity(proteins[i], proteins[j])
            if sim >= threshold:
                # Merge clusters
                ci, cj = clusters[i], clusters[j]
                for k in range(n):
                    if clusters[k] == cj:
                        clusters[k] = ci

    # Renumber clusters
    unique = sorted(set(clusters))
    mapping = {c: i for i, c in enumerate(unique)}
    clusters = [mapping[c] for c in clusters]
    return clusters

def main():
    data_dir = os.path.join(os.path.dirname(__file__), '..', 'data')

    # Read integrase DNA sequences
    int_fasta = os.path.join(data_dir, 'integrase_seqs.fasta')
    dna_seqs = read_fasta(int_fasta)

    # Translate to protein
    proteins = []
    names = []
    for name, dna in dna_seqs.items():
        prot = translate(dna)
        prot = prot.split('*')[0]  # stop at first stop codon
        if len(prot) >= 80:
            proteins.append(prot)
            names.append(name)

    print(f"Translated {len(proteins)} integrase proteins (≥80aa)")

    # Cluster
    clusters = cluster_proteins(proteins, names, threshold=0.30)
    n_clusters = max(clusters) + 1

    # Show clusters
    cluster_members = defaultdict(list)
    for i, c in enumerate(clusters):
        cluster_members[c].append(i)

    print(f"\nClusters (threshold 0.30): {n_clusters} families")
    for c in sorted(cluster_members.keys()):
        members = cluster_members[c]
        sizes = [len(proteins[i]) for i in members]
        print(f"  Family {c}: {len(members)} members, median size={statistics.median(sizes):.0f}aa")
        for i in members[:3]:
            desc = names[i][:60]
            print(f"    {desc} ({len(proteins[i])}aa)")
        if len(members) > 3:
            print(f"    ... and {len(members) - 3} more")

    # Build per-family profiles and test
    print(f"\n{'='*60}")
    print(f"  PER-FAMILY SCORING TEST")
    print(f"{'='*60}")

    # Collect all reference gene proteins for testing
    from integrase_profile import extract_marker_genes
    ref_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    family_profiles = []
    for c in sorted(cluster_members.keys()):
        members = cluster_members[c]
        if len(members) < 3:
            continue  # skip tiny families

        family_prots = [proteins[i] for i in members]
        profile = protein_3mer_profile(family_prots)
        family_profiles.append((c, len(members), profile))

    print(f"\n{len(family_profiles)} families with ≥3 members")

    # Test on E. coli: score all reference genes
    for genome in ["ecoli_k12", "salmonella_lt2"]:
        gff = f"{ref_dir}/{genome}.gff"
        fasta = f"{ref_dir}/{genome}.fasta"
        if not os.path.exists(gff): continue

        all_genes = extract_marker_genes(gff, fasta, [''])  # all CDS
        # Also get integrases specifically
        integrases = extract_marker_genes(gff, fasta, ['integrase'])

        print(f"\n--- {genome} ({len(all_genes)} genes, {len(integrases)} integrases) ---")

        # Score all genes with max across families
        int_set = set((g['start'], g['end']) for g in integrases)

        int_max_scores = []
        core_max_scores = []

        for g in all_genes:
            if g['length'] < 300: continue
            prot = translate(g['seq']).split('*')[0]
            if len(prot) < 50: continue

            max_score = max((score_protein(prot, prof) for _, _, prof in family_profiles), default=0)

            if (g['start'], g['end']) in int_set:
                int_max_scores.append(max_score)
            else:
                core_max_scores.append(max_score)

        if int_max_scores and core_max_scores:
            int_med = statistics.median(int_max_scores)
            core_med = statistics.median(core_max_scores)
            core_p95 = sorted(core_max_scores)[len(core_max_scores) * 95 // 100]
            core_p99 = sorted(core_max_scores)[len(core_max_scores) * 99 // 100]
            above_95 = sum(1 for s in int_max_scores if s > core_p95)
            above_99 = sum(1 for s in int_max_scores if s > core_p99)

            print(f"  Integrases max-family score: median={int_med:.4f}")
            print(f"  Core genes max-family score: median={core_med:.4f}")
            print(f"  Core 95th pct: {core_p95:.4f}, 99th: {core_p99:.4f}")
            print(f"  Integrases above core-95th: {above_95}/{len(int_max_scores)} ({above_95*100//len(int_max_scores)}%)")
            print(f"  Integrases above core-99th: {above_99}/{len(int_max_scores)} ({above_99*100//len(int_max_scores)}%)")

            # Threshold sweep
            for target_recall in [0.3, 0.5, 0.7]:
                target_tp = int(len(int_max_scores) * target_recall)
                all_sc = sorted([(s, 1) for s in int_max_scores] + [(s, 0) for s in core_max_scores], key=lambda x: -x[0])
                tp = fp = 0
                for score, label in all_sc:
                    if label == 1: tp += 1
                    else: fp += 1
                    if tp >= target_tp:
                        fpr = fp / len(core_max_scores)
                        print(f"  At {target_recall:.0%} recall: FPR={fpr:.2%} ({fp} FP)")
                        break

if __name__ == "__main__":
    main()
