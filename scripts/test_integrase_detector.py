#!/usr/bin/env python3
"""Test integrase detector: build hex profile from known integrases,
score all predicted genes, check separation."""

import sys, os, math, statistics

def read_fasta_seqs(path):
    seqs = {}
    name = ""
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name: seqs[name] = ''.join(parts)
                name = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip().upper())
    if name: seqs[name] = ''.join(parts)
    return seqs

def hex_enc(s):
    """Encode 6-mer to index 0-4095."""
    nt = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    v = 0
    for c in s:
        if c not in nt: return None
        v = v * 4 + nt[c]
    return v

def build_hex_profile(sequences):
    """Build log-odds hex profile from sequences vs uniform background."""
    counts = [0] * 4096
    total = 0
    for seq in sequences:
        for i in range(len(seq) - 5):
            idx = hex_enc(seq[i:i+6])
            if idx is not None:
                counts[idx] += 1
                total += 1

    # Background: uniform
    bg = 1.0 / 4096
    profile = [0.0] * 4096
    for i in range(4096):
        freq = (counts[i] + 0.5) / (total + 4096 * 0.5)
        profile[i] = math.log(freq / bg) if freq > 0 and bg > 0 else 0.0
    return profile

def score_sequence(seq, profile):
    """Score a DNA sequence against a hex profile."""
    total = 0.0
    n = 0
    for i in range(len(seq) - 5):
        idx = hex_enc(seq[i:i+6])
        if idx is not None:
            total += profile[idx]
            n += 1
    return total / n if n > 0 else 0.0

def extract_gene_seqs(gff_path, fasta_path):
    """Extract all CDS sequences with annotations."""
    seqs = read_fasta_seqs(fasta_path)
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS': continue
            if 'pseudo=true' in fields[8]: continue

            seqid = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            attrs = fields[8].lower()

            is_integrase = 'integrase' in attrs
            is_terminase = 'terminase' in attrs
            is_phage = any(kw in attrs for kw in ['phage', 'prophage', 'capsid', 'portal', 'tail', 'baseplate'])

            if seqid not in seqs: continue
            seq = seqs[seqid][start:end]
            if strand == '-':
                comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                seq = ''.join(comp.get(c, 'N') for c in reversed(seq))

            genes.append({
                'seq': seq, 'length': len(seq),
                'is_integrase': is_integrase,
                'is_terminase': is_terminase,
                'is_phage': is_phage,
            })
    return genes

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    # Read integrase sequences
    int_seqs = []
    int_fasta = os.path.join(os.path.dirname(__file__), '..', 'data', 'integrase_seqs.fasta')
    if os.path.exists(int_fasta):
        for name, seq in read_fasta_seqs(int_fasta).items():
            int_seqs.append(seq)

    print(f"Building profile from {len(int_seqs)} integrase sequences ({sum(len(s) for s in int_seqs)} bp)")

    # Build profile
    profile = build_hex_profile(int_seqs)

    # Test on each genome: score all CDS, check separation
    genomes = ["ecoli_k12", "salmonella_lt2", "bacillus_subtilis", "saureus_nctc8325"]

    for name in genomes:
        gff = f"{data_dir}/{name}.gff"
        fasta = f"{data_dir}/{name}.fasta"
        if not os.path.exists(gff): continue

        genes = extract_gene_seqs(gff, fasta)

        # Score all genes
        integrase_scores = []
        terminase_scores = []
        phage_scores = []
        core_scores = []

        for g in genes:
            if g['length'] < 200: continue
            score = score_sequence(g['seq'], profile)

            if g['is_integrase']:
                integrase_scores.append(score)
            elif g['is_terminase']:
                terminase_scores.append(score)
            elif g['is_phage']:
                phage_scores.append(score)
            else:
                core_scores.append(score)

        print(f"\n{'='*60}")
        print(f"  {name}")
        print(f"{'='*60}")

        if integrase_scores:
            print(f"  Integrases ({len(integrase_scores)}): median={statistics.median(integrase_scores):.3f}, "
                  f"min={min(integrase_scores):.3f}")
        if terminase_scores:
            print(f"  Terminases ({len(terminase_scores)}): median={statistics.median(terminase_scores):.3f}")
        if phage_scores:
            print(f"  Other phage ({len(phage_scores)}): median={statistics.median(phage_scores):.3f}")
        if core_scores:
            print(f"  Core genes ({len(core_scores)}): median={statistics.median(core_scores):.3f}")

            # Separation quality
            core_p95 = sorted(core_scores)[len(core_scores) * 95 // 100]
            core_p99 = sorted(core_scores)[len(core_scores) * 99 // 100]

            if integrase_scores:
                above_p95 = sum(1 for s in integrase_scores if s > core_p95)
                above_p99 = sum(1 for s in integrase_scores if s > core_p99)
                print(f"\n  Separation:")
                print(f"    Core 95th pct: {core_p95:.3f}")
                print(f"    Core 99th pct: {core_p99:.3f}")
                print(f"    Integrases above core-95th: {above_p95}/{len(integrase_scores)} ({above_p95*100//len(integrase_scores)}%)")
                print(f"    Integrases above core-99th: {above_p99}/{len(integrase_scores)} ({above_p99*100//len(integrase_scores)}%)")

                # At what threshold do we catch 50% of integrases with <5% FP?
                all_scores = [(s, 1) for s in integrase_scores] + [(s, 0) for s in core_scores]
                all_scores.sort(key=lambda x: -x[0])

                for target_recall in [0.3, 0.5, 0.7]:
                    target_tp = int(len(integrase_scores) * target_recall)
                    tp = 0
                    fp = 0
                    for score, label in all_scores:
                        if label == 1: tp += 1
                        else: fp += 1
                        if tp >= target_tp:
                            fpr = fp / len(core_scores)
                            print(f"    At {target_recall:.0%} recall: threshold={score:.3f}, FPR={fpr:.3%} ({fp} FP)")
                            break

if __name__ == "__main__":
    main()
