#!/usr/bin/env python3
"""Test att repeats on DETECTED candidate islands (not known).
Key question: do TP candidates have repeats and FP don't?"""

import os, re, subprocess

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

def overlap_frac(s1, e1, s2, e2):
    ov = max(0, min(e1, e2) - max(s1, s2))
    return ov / max(e1 - s1, 1)

def find_best_repeat(genome, left_pos, right_pos, flank=2000):
    """Find longest exact direct repeat between ±flank of left_pos and right_pos."""
    left_start = max(0, left_pos - flank)
    left_end = min(len(genome), left_pos + flank)
    right_start = max(0, right_pos - flank)
    right_end = min(len(genome), right_pos + flank)

    left_seq = genome[left_start:left_end]
    right_seq = genome[right_start:right_end]

    best_k = 0
    best_left_abs = 0
    best_right_abs = 0
    best_kmer = ""

    for k in range(15, 81):
        found = False
        left_kmers = {}
        for i in range(len(left_seq) - k + 1):
            kmer = left_seq[i:i+k]
            if 'N' not in kmer:
                left_kmers[kmer] = i

        for i in range(len(right_seq) - k + 1):
            kmer = right_seq[i:i+k]
            if kmer in left_kmers:
                best_k = k
                best_left_abs = left_start + left_kmers[kmer]
                best_right_abs = right_start + i
                best_kmer = kmer
                found = True
                break
        if not found and best_k > 0:
            break

    return best_k, best_left_abs, best_right_abs, best_kmer

def parse_detected_islands(stderr):
    islands = []
    for line in stderr.split('\n'):
        m = re.match(r'\s*Island:\s+(\d+)\.\.(\d+)\s+\(([\d.]+)kb,\s+(\d+)\s+genes\)\s+score=([\d.]+).*?(\w+)$', line)
        if m:
            islands.append({
                'start': int(m.group(1)),
                'end': int(m.group(2)),
                'span_kb': float(m.group(3)),
                'n_genes': int(m.group(4)),
                'score': float(m.group(5)),
                'type': m.group(6),
            })
    return islands

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"
    binary = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                          "target/release/prokrustes")

    for genome_name in ["ecoli_k12", "salmonella_lt2", "bacillus_subtilis", "saureus_nctc8325"]:
        gff = f"{data_dir}/{genome_name}.gff"
        fasta = f"{data_dir}/{genome_name}.fasta"
        ncrna = f"{data_dir}/{genome_name}_ncrna.gff"
        if not os.path.exists(fasta): continue

        seqs = read_fasta(fasta)
        genome = list(seqs.values())[0]
        known = find_known_phage_islands(gff)

        ncrna_flag = f"--ncrna {ncrna}" if os.path.exists(ncrna) else ""
        result = subprocess.run(f"{binary} {fasta} {ncrna_flag}",
                               shell=True, capture_output=True, text=True, timeout=120)
        detected = parse_detected_islands(result.stderr)

        if not detected: continue

        print(f"\n{'='*70}")
        print(f"  {genome_name}: {len(detected)} detected, {len(known)} known")
        print(f"{'='*70}")

        tp_with_repeat = 0
        tp_without = 0
        fp_with_repeat = 0
        fp_without = 0

        for det in detected:
            is_tp = any(overlap_frac(det['start'], det['end'], s, e) > 0.3 for s, e in known)
            label = "TP" if is_tp else "FP"

            # Search for att repeat
            rep_len, rep_left, rep_right, rep_seq = find_best_repeat(
                genome, det['start'], det['end'], flank=2000)

            has_repeat = rep_len >= 15

            if is_tp:
                if has_repeat: tp_with_repeat += 1
                else: tp_without += 1
            else:
                if has_repeat: fp_with_repeat += 1
                else: fp_without += 1

            if has_repeat:
                print(f"  {det['start']:>8}..{det['end']:>8} ({det['span_kb']:.0f}kb) {label} "
                      f"REPEAT {rep_len}bp dist_L={rep_left-det['start']}bp dist_R={rep_right-det['end']}bp")
            else:
                print(f"  {det['start']:>8}..{det['end']:>8} ({det['span_kb']:.0f}kb) {label} no repeat ≥15bp")

        print(f"\n  Summary:")
        print(f"    TP with repeat:    {tp_with_repeat}")
        print(f"    TP without repeat: {tp_without}")
        print(f"    FP with repeat:    {fp_with_repeat}")
        print(f"    FP without repeat: {fp_without}")
        if tp_with_repeat + fp_with_repeat > 0:
            prec_with = tp_with_repeat / (tp_with_repeat + fp_with_repeat)
            print(f"    Precision WITH repeat:    {prec_with:.2f}")
        if tp_without + fp_without > 0:
            prec_without = tp_without / (tp_without + fp_without)
            print(f"    Precision WITHOUT repeat: {prec_without:.2f}")

if __name__ == "__main__":
    main()
