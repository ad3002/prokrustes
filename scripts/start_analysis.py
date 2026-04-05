#!/usr/bin/env python3
"""Analyze how many errors are start-selection errors (same stop, wrong start)."""
import sys
from evaluate import parse_prediction_gff, parse_reference_gff

def analyze_starts(pred_path, ref_path):
    predicted = parse_prediction_gff(pred_path)
    known = parse_reference_gff(ref_path)

    # For each prediction, find the known gene with the closest stop codon
    start_errors = 0
    exact_matches = 0
    close_matches = 0  # within 30bp
    total_matched = 0

    for pred in predicted:
        # Find known gene on same strand with same or very close end
        best_end_diff = float('inf')
        best_known = None
        for gene in known:
            if pred["strand"] != gene["strand"]:
                continue
            end_diff = abs(pred["end"] - gene["end"])
            if end_diff < best_end_diff:
                best_end_diff = end_diff
                best_known = gene

        if best_known is None:
            continue

        if best_end_diff <= 3:  # Same stop group (within 3bp)
            total_matched += 1
            start_diff = abs(pred["start"] - best_known["start"])
            if start_diff <= 3:
                exact_matches += 1
            elif start_diff <= 30:
                close_matches += 1
            else:
                start_errors += 1

    # FN: known genes with no matching prediction (same stop)
    fn_no_stop_match = 0
    fn_wrong_start = 0
    for gene in known:
        best_end_diff = float('inf')
        best_pred = None
        for pred in predicted:
            if pred["strand"] != gene["strand"]:
                continue
            end_diff = abs(pred["end"] - gene["end"])
            if end_diff < best_end_diff:
                best_end_diff = end_diff
                best_pred = pred

        if best_pred is None or best_end_diff > 3:
            fn_no_stop_match += 1
        else:
            start_diff = abs(best_pred["start"] - gene["start"])
            if start_diff > 30:
                fn_wrong_start += 1

    print(f"  Same-stop matches: {total_matched}")
    print(f"    Exact start (±3bp): {exact_matches} ({exact_matches*100//max(total_matched,1)}%)")
    print(f"    Close start (±30bp): {close_matches}")
    print(f"    Wrong start (>30bp): {start_errors}")
    print(f"  FN: no stop match={fn_no_stop_match}, wrong start={fn_wrong_start}")

if __name__ == "__main__":
    analyze_starts(sys.argv[1], sys.argv[2])
