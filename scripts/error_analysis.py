#!/usr/bin/env python3
"""Detailed error analysis: classify FP/FN by type."""
import sys
from evaluate import parse_prediction_gff, parse_reference_gff

def analyze(pred_path, ref_path):
    predicted = parse_prediction_gff(pred_path)
    known = parse_reference_gff(ref_path)

    # Match predictions to reference
    matched_ref = [None] * len(known)
    matched_pred = [None] * len(predicted)

    for pi, pred in enumerate(predicted):
        best_ov, best_idx = 0.0, None
        for ki, gene in enumerate(known):
            if matched_ref[ki] is not None:
                continue
            if pred["strand"] != gene["strand"]:
                continue
            ov_s = max(pred["start"], gene["start"])
            ov_e = min(pred["end"], gene["end"])
            if ov_s >= ov_e:
                continue
            ov_len = ov_e - ov_s
            pred_len = max(pred["end"] - pred["start"], 1)
            gene_len = max(gene["end"] - gene["start"], 1)
            recip = min(ov_len / pred_len, ov_len / gene_len)
            if recip > best_ov:
                best_ov, best_idx = recip, ki
        if best_ov >= 0.8 and best_idx is not None:
            matched_ref[best_idx] = pi
            matched_pred[pi] = best_idx

    tp = sum(1 for m in matched_ref if m is not None)

    # Classify FP
    fp_short = 0  # <300bp
    fp_medium = 0  # 300-600bp
    fp_long = 0  # >600bp
    fp_overlap = 0  # overlaps a known gene (wrong boundaries)
    fp_no_overlap = 0  # no overlap at all

    for pi, pred in enumerate(predicted):
        if matched_pred[pi] is not None:
            continue
        plen = pred["end"] - pred["start"]

        # Check if it partially overlaps any known gene
        has_partial = False
        for ki, gene in enumerate(known):
            if pred["strand"] != gene["strand"]:
                continue
            ov_s = max(pred["start"], gene["start"])
            ov_e = min(pred["end"], gene["end"])
            if ov_s < ov_e:
                has_partial = True
                break

        if has_partial:
            fp_overlap += 1
        else:
            fp_no_overlap += 1

        if plen < 300:
            fp_short += 1
        elif plen < 600:
            fp_medium += 1
        else:
            fp_long += 1

    # Classify FN
    fn_short = 0  # <300bp
    fn_medium = 0  # 300-600bp
    fn_long = 0  # >600bp
    fn_partial = 0  # partially overlapped by a prediction
    fn_invisible = 0  # no prediction overlaps at all

    for ki, gene in enumerate(known):
        if matched_ref[ki] is not None:
            continue
        glen = gene["end"] - gene["start"]

        has_partial = False
        for pi, pred in enumerate(predicted):
            if pred["strand"] != gene["strand"]:
                continue
            ov_s = max(pred["start"], gene["start"])
            ov_e = min(pred["end"], gene["end"])
            if ov_s < ov_e:
                has_partial = True
                break

        if has_partial:
            fn_partial += 1
        else:
            fn_invisible += 1

        if glen < 300:
            fn_short += 1
        elif glen < 600:
            fn_medium += 1
        else:
            fn_long += 1

    total_fp = fp_short + fp_medium + fp_long
    total_fn = fn_short + fn_medium + fn_long

    print(f"  TP={tp}  FP={total_fp}  FN={total_fn}")
    print(f"  FP breakdown: short(<300)={fp_short} med(300-600)={fp_medium} long(>600)={fp_long}")
    print(f"  FP type: overlap_wrong_boundary={fp_overlap} no_overlap={fp_no_overlap}")
    print(f"  FN breakdown: short(<300)={fn_short} med(300-600)={fn_medium} long(>600)={fn_long}")
    print(f"  FN type: partial_overlap={fn_partial} invisible={fn_invisible}")

if __name__ == "__main__":
    analyze(sys.argv[1], sys.argv[2])
