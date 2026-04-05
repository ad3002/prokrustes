#!/usr/bin/env python3
"""Analyze FP: are they short genes that shadow real ones?"""
import sys
from evaluate import parse_prediction_gff, parse_reference_gff

def analyze_fp(pred_path, ref_path):
    predicted = parse_prediction_gff(pred_path)
    known = parse_reference_gff(ref_path)

    # Match
    matched_pred = set()
    matched_ref = set()
    for pi, pred in enumerate(predicted):
        best_ov, best_idx = 0.0, None
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
            if recip > best_ov:
                best_ov, best_idx = recip, ki
        if best_ov >= 0.8 and best_idx is not None:
            matched_pred.add(pi)
            matched_ref.add(best_idx)

    # Analyze FP
    fp_near_known = 0
    fp_shadows_known = 0  # FP on opposite strand of a real gene
    fp_isolated = 0
    fp_by_len = {"<150": 0, "150-300": 0, "300-600": 0, ">600": 0}

    for pi, pred in enumerate(predicted):
        if pi in matched_pred: continue
        plen = pred["end"] - pred["start"]

        if plen < 150: fp_by_len["<150"] += 1
        elif plen < 300: fp_by_len["150-300"] += 1
        elif plen < 600: fp_by_len["300-600"] += 1
        else: fp_by_len[">600"] += 1

        # Check if it's near/overlapping a known gene
        near = False
        shadow = False
        for ki, gene in enumerate(known):
            ov_s = max(pred["start"], gene["start"])
            ov_e = min(pred["end"], gene["end"])
            if ov_e > ov_s:
                near = True
                if pred["strand"] != gene["strand"]:
                    shadow = True

        if shadow: fp_shadows_known += 1
        elif near: fp_near_known += 1
        else: fp_isolated += 1

    # Analyze FN
    fn_by_len = {"<150": 0, "150-300": 0, "300-600": 0, ">600": 0}
    fn_has_partial = 0

    for ki, gene in enumerate(known):
        if ki in matched_ref: continue
        glen = gene["end"] - gene["start"]

        if glen < 150: fn_by_len["<150"] += 1
        elif glen < 300: fn_by_len["150-300"] += 1
        elif glen < 600: fn_by_len["300-600"] += 1
        else: fn_by_len[">600"] += 1

        # Check partial overlap
        for pi, pred in enumerate(predicted):
            if pred["strand"] != gene["strand"]: continue
            ov_s = max(pred["start"], gene["start"])
            ov_e = min(pred["end"], gene["end"])
            if ov_e > ov_s:
                fn_has_partial += 1
                break

    total_fp = len(predicted) - len(matched_pred)
    total_fn = len(known) - len(matched_ref)

    print(f"  FP={total_fp}: near_known={fp_near_known} shadow_opp={fp_shadows_known} isolated={fp_isolated}")
    print(f"  FP by length: {fp_by_len}")
    print(f"  FN={total_fn}: has_partial_overlap={fn_has_partial} invisible={total_fn-fn_has_partial}")
    print(f"  FN by length: {fn_by_len}")

if __name__ == "__main__":
    analyze_fp(sys.argv[1], sys.argv[2])
