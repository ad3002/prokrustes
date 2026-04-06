#!/usr/bin/env python3
"""Calibrate prophage scoring law from labeled data.
Collects segment features from all genomes, labels TP/FP,
fits optimal weights and threshold."""

import sys, os, re, subprocess, json

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

def parse_island_features(stderr):
    """Parse Island: lines with features from Prokrustes stderr."""
    islands = []
    for line in stderr.split('\n'):
        m = re.match(r'\s*Island:\s+(\d+)\.\.(\d+)\s+\(([\d.]+)kb,\s+(\d+)\s+genes\)\s+score=([\d.]+)\s+\[F=([\d.]+)\s+D=([\d.]+)\s+frag=([\d.]+)\s+maxRun=(\d+)\]\s+(\w+)', line)
        if m:
            islands.append({
                'start': int(m.group(1)),
                'end': int(m.group(2)),
                'span_kb': float(m.group(3)),
                'n_genes': int(m.group(4)),
                'score': float(m.group(5)),
                'F': float(m.group(6)),  # lowHexFrac evidence
                'D': float(m.group(7)),  # densityExcess
                'frag': float(m.group(8)),  # fragmentation
                'maxRun': int(m.group(9)),
                'type': m.group(10),
            })
    return islands

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"
    binary = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                          "target/release/prokrustes")

    genomes = ["ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
               "saureus_nctc8325", "bacteroides_fragilis", "borrelia_burgdorferi"]

    all_tp_features = []
    all_fp_features = []

    for name in genomes:
        fasta = f"{data_dir}/{name}.fasta"
        ref = f"{data_dir}/{name}.gff"
        ncrna = f"{data_dir}/{name}_ncrna.gff"
        if not os.path.exists(fasta) or not os.path.exists(ref): continue

        known = find_known_phage_islands(ref)
        if not known: continue

        ncrna_flag = f"--ncrna {ncrna}" if os.path.exists(ncrna) else ""
        result = subprocess.run(f"{binary} {fasta} {ncrna_flag}",
                               shell=True, capture_output=True, text=True, timeout=120)

        detected = parse_island_features(result.stderr)

        for det in detected:
            is_tp = any(overlap_frac(det['start'], det['end'], s, e) > 0.3
                       for s, e in known)
            entry = {
                'genome': name,
                'F': det['F'],
                'D': det['D'],
                'frag': det['frag'],
                'maxRun': det['maxRun'],
                'score': det['score'],
                'span_kb': det['span_kb'],
                'n_genes': det['n_genes'],
                'label': 1 if is_tp else 0,
            }
            if is_tp:
                all_tp_features.append(entry)
            else:
                all_fp_features.append(entry)

    print(f"Collected: {len(all_tp_features)} TP segments, {len(all_fp_features)} FP segments")
    print()

    # Show feature distributions
    for feat in ['F', 'D', 'frag', 'maxRun', 'score', 'span_kb', 'n_genes']:
        tp_vals = sorted([e[feat] for e in all_tp_features])
        fp_vals = sorted([e[feat] for e in all_fp_features])
        if not tp_vals or not fp_vals: continue

        import statistics
        tp_med = statistics.median(tp_vals)
        fp_med = statistics.median(fp_vals)
        tp_p25 = tp_vals[len(tp_vals)//4]
        fp_p75 = fp_vals[len(fp_vals)*3//4]

        print(f"{feat:>10}: TP median={tp_med:.3f} (25th={tp_p25:.3f}), FP median={fp_med:.3f} (75th={fp_p75:.3f})")

    print()

    # Try different thresholds on the current score
    print("=== Current score threshold sweep ===")
    all_entries = all_tp_features + all_fp_features
    for thresh in [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 1.7, 2.0]:
        tp = sum(1 for e in all_tp_features if e['score'] >= thresh)
        fp = sum(1 for e in all_fp_features if e['score'] >= thresh)
        fn = len(all_tp_features) - tp
        prec = tp / max(tp + fp, 1)
        rec = tp / max(tp + fn, 1)
        f1 = 2 * prec * rec / max(prec + rec, 0.001)
        print(f"  thresh={thresh:.1f}: TP={tp:3d} FP={fp:3d} FN={fn:2d} prec={prec:.2f} rec={rec:.2f} F1={f1:.2f}")

    print()

    # Grid search over weights
    print("=== Weight optimization (grid search) ===")
    best_f1 = 0
    best_params = None

    for w1_10 in range(5, 25, 3):        # lowHexFrac
        for w3_10 in range(0, 10, 2):     # density
            for w4_10 in range(2, 15, 2): # frag penalty
                for w5_10 in range(2, 12, 2):  # size (span_kb)
                    for w6_10 in range(2, 12, 2):  # n_genes
                        w1 = w1_10 / 10.0
                        w3 = w3_10 / 10.0
                        w4 = w4_10 / 10.0
                        w5 = w5_10 / 10.0
                        w6 = w6_10 / 10.0

                        scores = []
                        for e in all_entries:
                            f = min(1.0, e['F'] / 0.5)
                            d = min(1.0, e['D'])
                            sz = min(1.5, e['span_kb'] / 20.0)
                            ng = min(1.5, e['n_genes'] / 20.0)
                            s = w1*f + w3*d + w5*sz + w6*ng - w4*e['frag']
                            scores.append((s, e['label']))

                    # Find best threshold for this weight combo
                    for thresh_10 in range(3, 25, 1):
                        thresh = thresh_10 / 10.0
                        tp = sum(1 for s, l in scores if s >= thresh and l == 1)
                        fp = sum(1 for s, l in scores if s >= thresh and l == 0)
                        fn = sum(1 for s, l in scores if s < thresh and l == 1)
                        prec = tp / max(tp + fp, 1)
                        rec = tp / max(tp + fn, 1)
                        f1 = 2 * prec * rec / max(prec + rec, 0.001)

                        if f1 > best_f1:
                            best_f1 = f1
                            best_params = (w1, w3, w4, w5, w6, thresh, tp, fp, fn, prec, rec)

    if best_params:
        w1, w3, w4, w5, w6, thresh, tp, fp, fn, prec, rec = best_params
        print(f"  Best: w_F={w1:.1f} w_D={w3:.1f} w_frag={w4:.1f} w_size={w5:.1f} w_ngenes={w6:.1f} thresh={thresh:.1f}")
        print(f"  TP={tp} FP={fp} FN={fn} prec={prec:.2f} rec={rec:.2f} F1={best_f1:.2f}")

if __name__ == "__main__":
    main()
