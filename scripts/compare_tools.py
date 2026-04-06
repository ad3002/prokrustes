#!/usr/bin/env python3
"""Compare Prokrustes vs Prodigal predictions directly.
Find genes that one tool finds and the other misses."""
import sys
from evaluate import parse_prediction_gff, parse_reference_gff

def match(preds, refs):
    """Match predictions to reference. Returns (tp_pred_idx, tp_ref_idx, fp_pred_idx, fn_ref_idx)."""
    tp_p, tp_r = set(), set()
    for pi, pred in enumerate(preds):
        for ri, ref in enumerate(refs):
            if ri in tp_r: continue
            if pred["strand"] != ref["strand"]: continue
            ov_s = max(pred["start"], ref["start"])
            ov_e = min(pred["end"], ref["end"])
            if ov_s >= ov_e: continue
            ov = ov_e - ov_s
            pl = max(pred["end"] - pred["start"], 1)
            rl = max(ref["end"] - ref["start"], 1)
            if min(ov/pl, ov/rl) >= 0.8:
                tp_p.add(pi)
                tp_r.add(ri)
                break
    fp = set(range(len(preds))) - tp_p
    fn = set(range(len(refs))) - tp_r
    return tp_p, tp_r, fp, fn

def main():
    pk_path, pd_path, ref_path = sys.argv[1], sys.argv[2], sys.argv[3]

    pk = parse_prediction_gff(pk_path)
    pd = parse_prediction_gff(pd_path)
    ref = parse_reference_gff(ref_path)

    pk_tp, pk_tp_r, pk_fp, pk_fn = match(pk, ref)
    pd_tp, pd_tp_r, pd_fp, pd_fn = match(pd, ref)

    # Genes Prodigal finds that we miss
    pd_finds_we_miss = pd_tp_r - pk_tp_r
    # Genes we find that Prodigal misses
    we_find_pd_misses = pk_tp_r - pd_tp_r

    print(f"Prokrustes: TP={len(pk_tp)} FP={len(pk_fp)} FN={len(pk_fn)}")
    print(f"Prodigal:   TP={len(pd_tp)} FP={len(pd_fp)} FN={len(pd_fn)}")
    print(f"Prodigal finds, we miss: {len(pd_finds_we_miss)}")
    print(f"We find, Prodigal misses: {len(we_find_pd_misses)}")
    print()

    # Characterize the genes Prodigal finds that we miss
    if pd_finds_we_miss:
        lens = [ref[ri]["end"] - ref[ri]["start"] for ri in pd_finds_we_miss]
        lens.sort()
        print(f"Genes Prodigal finds that we miss ({len(lens)}):")
        print(f"  Length distribution: min={min(lens)} median={lens[len(lens)//2]} max={max(lens)}")
        short = sum(1 for l in lens if l < 300)
        med = sum(1 for l in lens if 300 <= l < 600)
        lng = sum(1 for l in lens if l >= 600)
        print(f"  <300bp: {short}, 300-600bp: {med}, >600bp: {lng}")

if __name__ == "__main__":
    main()
