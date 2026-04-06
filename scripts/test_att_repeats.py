#!/usr/bin/env python3
"""Test: do known phage islands have direct repeats (attL/attR) at boundaries?
For each known island, extract ±5kb flanks and find shared exact k-mers."""

import os, statistics

def read_fasta(path):
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

def find_shared_kmers(left_seq, right_seq, k):
    """Find k-mers shared between left and right boundary sequences."""
    if len(left_seq) < k or len(right_seq) < k:
        return []

    left_kmers = {}
    for i in range(len(left_seq) - k + 1):
        kmer = left_seq[i:i+k]
        if 'N' not in kmer:
            if kmer not in left_kmers:
                left_kmers[kmer] = []
            left_kmers[kmer].append(i)

    matches = []
    for i in range(len(right_seq) - k + 1):
        kmer = right_seq[i:i+k]
        if kmer in left_kmers:
            for left_pos in left_kmers[kmer]:
                matches.append((kmer, left_pos, i))

    return matches

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"
    flank = 5000  # 5kb flanks

    for genome_name in ["ecoli_k12", "salmonella_lt2", "bacillus_subtilis", "saureus_nctc8325"]:
        gff = f"{data_dir}/{genome_name}.gff"
        fasta = f"{data_dir}/{genome_name}.fasta"
        if not os.path.exists(gff): continue

        seqs = read_fasta(fasta)
        genome = list(seqs.values())[0]  # first/main chromosome
        islands = find_known_phage_islands(gff)
        if not islands: continue

        print(f"\n{'='*70}")
        print(f"  {genome_name} ({len(islands)} known islands)")
        print(f"{'='*70}")

        for s, e in islands:
            size_kb = (e - s) / 1000
            if size_kb < 5: continue  # skip tiny islands

            # Extract boundary regions
            left_start = max(0, s - flank)
            left_end = min(len(genome), s + flank)
            right_start = max(0, e - flank)
            right_end = min(len(genome), e + flank)

            left_seq = genome[left_start:left_end]
            right_seq = genome[right_start:right_end]

            # Search for shared direct repeats at various lengths
            best_k = 0
            best_match = None
            best_count = 0

            for k in range(10, 81):
                matches = find_shared_kmers(left_seq, right_seq, k)
                if matches:
                    best_k = k
                    best_match = matches[0]
                    best_count = len(matches)

            if best_match:
                kmer, lpos, rpos = best_match
                # Absolute positions
                abs_left = left_start + lpos
                abs_right = right_start + rpos
                print(f"  Island {s}..{e} ({size_kb:.0f}kb): "
                      f"REPEAT found! {best_k}bp at left={abs_left} right={abs_right}")
                print(f"    Sequence: {kmer[:40]}{'...' if len(kmer) > 40 else ''}")
                print(f"    Distance from boundary: left={abs_left-s}bp, right={abs_right-e}bp")
            else:
                # Try shorter k-mers
                for k in [8, 9]:
                    matches = find_shared_kmers(left_seq, right_seq, k)
                    if len(matches) > 20:  # many short matches = noise
                        continue
                    if matches:
                        kmer, lpos, rpos = matches[0]
                        abs_left = left_start + lpos
                        abs_right = right_start + rpos
                        print(f"  Island {s}..{e} ({size_kb:.0f}kb): "
                              f"weak repeat {k}bp ({len(matches)} hits)")
                        break
                else:
                    print(f"  Island {s}..{e} ({size_kb:.0f}kb): NO repeat found (≥8bp)")

        # Also test: do random non-island positions have repeats?
        print(f"\n  Control: random positions")
        import random
        random.seed(42)
        n_with_repeat = 0
        n_tested = 20
        for _ in range(n_tested):
            pos = random.randint(flank, len(genome) - flank - 20000)
            fake_end = pos + 20000
            left_seq = genome[pos-flank:pos+flank]
            right_seq = genome[fake_end-flank:fake_end+flank]
            for k in range(10, 81):
                matches = find_shared_kmers(left_seq, right_seq, k)
                if matches:
                    n_with_repeat += 1
                    break
        print(f"    {n_with_repeat}/{n_tested} random 20kb segments have ≥10bp shared repeat")

if __name__ == "__main__":
    main()
