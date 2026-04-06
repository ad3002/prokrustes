#!/usr/bin/env python3
"""Analyze start-site errors: how many genes have correct stop but wrong start?
What features distinguish correct vs incorrect starts?"""

import os, sys, statistics
sys.path.insert(0, os.path.dirname(__file__))
from evaluate import parse_prediction_gff, parse_reference_gff

def analyze_starts(pred_path, ref_path):
    preds = parse_prediction_gff(pred_path)
    refs = parse_reference_gff(ref_path)

    # Match by stop codon (end position ±3bp, same strand)
    exact_start = 0
    close_start = 0     # ±30bp
    wrong_start = 0      # >30bp
    too_long = 0         # our start upstream of reference (gene too long)
    too_short = 0        # our start downstream of reference (gene too short)
    no_stop_match = 0

    wrong_start_deltas = []  # signed: positive = too long, negative = too short

    for pred in preds:
        best_ref = None
        best_end_diff = float('inf')
        for ref in refs:
            if pred["strand"] != ref["strand"]: continue
            end_diff = abs(pred["end"] - ref["end"])
            if end_diff < best_end_diff:
                best_end_diff = end_diff
                best_ref = ref

        if best_ref is None or best_end_diff > 3:
            no_stop_match += 1
            continue

        start_diff = pred["start"] - best_ref["start"]  # positive = our start is downstream (too short)

        if abs(start_diff) <= 3:
            exact_start += 1
        elif abs(start_diff) <= 30:
            close_start += 1
        else:
            wrong_start += 1
            wrong_start_deltas.append(start_diff)
            if start_diff > 0:
                too_short += 1
            else:
                too_long += 1

    return {
        "exact": exact_start,
        "close": close_start,
        "wrong": wrong_start,
        "too_long": too_long,
        "too_short": too_short,
        "no_stop": no_stop_match,
        "deltas": wrong_start_deltas,
    }

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    total_exact = 0
    total_wrong = 0
    total_long = 0
    total_short = 0
    all_deltas = []

    for name in ["ecoli_k12", "salmonella_lt2", "bacillus_subtilis", "saureus_nctc8325",
                  "mycoplasma_genitalium", "borrelia_burgdorferi", "bacteroides_fragilis",
                  "synechocystis_pcc6803", "pseudomonas_aeruginosa", "mycobacterium_tuberculosis"]:
        pred_path = f"/tmp/{name}_test.gff"
        ref_path = f"{data_dir}/{name}.gff"
        if not os.path.exists(pred_path): continue

        r = analyze_starts(pred_path, ref_path)

        total_matched = r["exact"] + r["close"] + r["wrong"]
        accuracy = r["exact"] / max(total_matched, 1)

        print(f"{name}:")
        print(f"  Same-stop matches: {total_matched}")
        print(f"  Exact start (±3bp): {r['exact']} ({accuracy:.1%})")
        print(f"  Close (±30bp): {r['close']}")
        print(f"  Wrong (>30bp): {r['wrong']} (too_long={r['too_long']}, too_short={r['too_short']})")
        if r['deltas']:
            abs_deltas = [abs(d) for d in r['deltas']]
            print(f"  Wrong start offset: median={statistics.median(abs_deltas):.0f}bp, max={max(abs_deltas)}bp")
        print(f"  No stop match: {r['no_stop']}")
        print()

        total_exact += r["exact"]
        total_wrong += r["wrong"]
        total_long += r["too_long"]
        total_short += r["too_short"]
        all_deltas.extend(r["deltas"])

    print(f"{'='*60}")
    print(f"TOTAL across all genomes:")
    print(f"  Exact starts: {total_exact}")
    print(f"  Wrong starts: {total_wrong} (too_long={total_long}, too_short={total_short})")
    if all_deltas:
        abs_all = [abs(d) for d in all_deltas]
        print(f"  Median offset: {statistics.median(abs_all):.0f}bp")
        for bucket in [30, 60, 100, 200, 500, 1000]:
            n = sum(1 for d in abs_all if d <= bucket)
            print(f"    ≤{bucket}bp: {n} ({n*100//len(abs_all)}%)")

    print(f"\n  Theoretical F1 gain if ALL wrong starts fixed:")
    print(f"    Each wrong start is both 1 FP + 1 FN")
    print(f"    Fixing {total_wrong} wrong starts = -{total_wrong} FP, -{total_wrong} FN")

if __name__ == "__main__":
    main()
