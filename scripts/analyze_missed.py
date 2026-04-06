#!/usr/bin/env python3
"""Find genes Prodigal finds that we miss, check what our pipeline scored them."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from evaluate import parse_prediction_gff, parse_reference_gff

def match_to_ref(preds, refs):
    tp_r = set()
    for pi, p in enumerate(preds):
        for ri, r in enumerate(refs):
            if ri in tp_r: continue
            if p["strand"] != r["strand"]: continue
            ov_s = max(p["start"], r["start"])
            ov_e = min(p["end"], r["end"])
            if ov_s >= ov_e: continue
            ov = ov_e - ov_s
            if min(ov/max(p["end"]-p["start"],1), ov/max(r["end"]-r["start"],1)) >= 0.8:
                tp_r.add(ri)
                break
    return tp_r

def main():
    pk_path, pd_path, ref_path = sys.argv[1], sys.argv[2], sys.argv[3]
    all_orfs_path = sys.argv[4] if len(sys.argv) > 4 else None

    pk = parse_prediction_gff(pk_path)
    pd = parse_prediction_gff(pd_path)
    refs = parse_reference_gff(ref_path)

    pk_matched = match_to_ref(pk, refs)
    pd_matched = match_to_ref(pd, refs)

    # Genes Prodigal finds that we miss
    pd_only = pd_matched - pk_matched

    print(f"Prodigal finds {len(pd_only)} genes that we miss:")
    print()

    # Parse all-orfs dump if provided
    all_orfs = {}
    if all_orfs_path and os.path.exists(all_orfs_path):
        with open(all_orfs_path) as f:
            header = f.readline().strip().split('\t')
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 17:
                    continue
                start = int(parts[0])
                end = int(parts[1])
                strand = parts[2]
                length = int(parts[3])
                hex_avg = float(parts[6])
                rbs = float(parts[10])
                score = float(parts[16])
                key = (start, end, strand)
                all_orfs[key] = {
                    "length": length, "hex_avg": hex_avg,
                    "rbs": rbs, "score": score,
                }

    for ri in sorted(pd_only):
        ref = refs[ri]
        rlen = ref["end"] - ref["start"]
        # Check if we had this ORF in all_orfs
        orf_info = ""
        if all_orfs:
            key = (ref["start"], ref["end"], ref["strand"])
            if key in all_orfs:
                o = all_orfs[key]
                orf_info = f" [found in ORFs: hex={o['hex_avg']:.3f} rbs={o['rbs']:.2f} score={o['score']:.3f}]"
            else:
                # Try nearby
                for (s, e, st), o in all_orfs.items():
                    if st != ref["strand"]: continue
                    if abs(e - ref["end"]) <= 3 and abs(s - ref["start"]) <= 100:
                        orf_info = f" [near-match {s}..{e}: hex={o['hex_avg']:.3f} rbs={o['rbs']:.2f} score={o['score']:.3f}]"
                        break
                if not orf_info:
                    orf_info = " [NOT FOUND in ORFs]"

        print(f"  {ref['start']}..{ref['end']} ({ref['strand']}) {rlen}bp{orf_info}")

if __name__ == "__main__":
    main()
