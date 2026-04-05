#!/usr/bin/env python3
"""Learn what distinguishes TP from FP predictions across all genomes.
Uses prediction GFF (with feature attributes) and reference GFF."""

import sys, os, re
from evaluate import parse_reference_gff

def parse_prokrustes_gff(path):
    """Parse Prokrustes GFF with attributes."""
    genes = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS': continue
            attrs = {}
            for part in fields[8].split(';'):
                if '=' in part:
                    k, v = part.split('=', 1)
                    attrs[k] = v
            genes.append({
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "score": float(fields[5]),
                "length": int(fields[4]) - int(fields[3]) + 1,
                "start_type": attrs.get("start_type", "?"),
                "rbs": float(attrs.get("rbs_score", "0")),
                "hex": float(attrs.get("hex_score", "0")),
            })
    return genes

def match_predictions(predicted, known):
    """Return sets of TP and FP indices."""
    tp_set = set()
    matched_ref = set()
    for pi, pred in enumerate(predicted):
        for ki, gene in enumerate(known):
            if ki in matched_ref: continue
            if pred["strand"] != gene["strand"]: continue
            ov_s = max(pred["start"], gene["start"])
            ov_e = min(pred["end"], gene["end"])
            if ov_s >= ov_e: continue
            ov_len = ov_e - ov_s
            pred_len = max(pred["end"] - pred["start"], 1)
            gene_len = max(gene["end"] - gene["start"], 1)
            recip = min(ov_len / pred_len, ov_len / gene_len)
            if recip >= 0.8:
                tp_set.add(pi)
                matched_ref.add(ki)
                break
    return tp_set

def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "/Users/akomissarov/Dropbox/workspace/new/renorm/data"
    pred_dir = sys.argv[2] if len(sys.argv) > 2 else "/tmp"

    genomes = [
        "ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
        "pseudomonas_aeruginosa", "saureus_nctc8325",
        "mycobacterium_tuberculosis", "mycoplasma_genitalium",
        "synechocystis_pcc6803", "bacteroides_fragilis",
        "borrelia_burgdorferi"
    ]

    all_tp = []
    all_fp = []

    for name in genomes:
        pred_path = os.path.join(pred_dir, f"{name}_test.gff")
        ref_path = os.path.join(data_dir, f"{name}.gff")
        if not os.path.exists(pred_path) or not os.path.exists(ref_path):
            continue

        predicted = parse_prokrustes_gff(pred_path)
        known = parse_reference_gff(ref_path)
        tp_set = match_predictions(predicted, known)

        for pi, pred in enumerate(predicted):
            entry = {
                "genome": name,
                "length": pred["length"],
                "score": pred["score"],
                "rbs": pred["rbs"],
                "hex": pred["hex"],
                "start_type": pred["start_type"],
            }
            if pi in tp_set:
                all_tp.append(entry)
            else:
                all_fp.append(entry)

    print(f"Total: {len(all_tp)} TP, {len(all_fp)} FP")
    print()

    # Compare distributions
    for feature in ["length", "score", "rbs", "hex"]:
        tp_vals = [e[feature] for e in all_tp]
        fp_vals = [e[feature] for e in all_fp]
        if not tp_vals or not fp_vals: continue

        import statistics
        tp_med = statistics.median(tp_vals)
        fp_med = statistics.median(fp_vals)
        tp_p10 = sorted(tp_vals)[len(tp_vals)//10]
        fp_p90 = sorted(fp_vals)[len(fp_vals)*9//10]

        print(f"{feature}:")
        print(f"  TP: median={tp_med:.3f}, 10th pct={tp_p10:.3f}")
        print(f"  FP: median={fp_med:.3f}, 90th pct={fp_p90:.3f}")

        # Threshold that captures 95% of TP
        tp_p5 = sorted(tp_vals)[len(tp_vals)//20]
        print(f"  TP 5th pct (threshold)={tp_p5:.3f}")
        fp_below = sum(1 for v in fp_vals if v < tp_p5)
        print(f"  FP below TP-5th: {fp_below}/{len(fp_vals)} ({fp_below*100//len(fp_vals)}%)")
        print()

    # Length bucketed analysis
    for bucket, lo, hi in [("<150", 0, 150), ("150-300", 150, 300), ("300-600", 300, 600), (">600", 600, 99999)]:
        tp_n = sum(1 for e in all_tp if lo <= e["length"] < hi)
        fp_n = sum(1 for e in all_fp if lo <= e["length"] < hi)
        total = tp_n + fp_n
        ratio = tp_n / max(total, 1)
        print(f"  {bucket}: {tp_n} TP, {fp_n} FP, precision={ratio:.2f}")

if __name__ == "__main__":
    main()
