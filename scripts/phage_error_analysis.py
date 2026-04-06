#!/usr/bin/env python3
"""Analyze TP/FP/FN specifically in prophage regions vs core genome."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from evaluate import parse_prediction_gff, parse_reference_gff

def find_phage_regions(gff_path):
    """Find prophage island boundaries from NCBI GFF phage-related annotations."""
    phage_genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9: continue
            if fields[2] != 'CDS': continue
            attrs = fields[8].lower()
            if any(kw in attrs for kw in ['phage', 'prophage', 'capsid', 'terminase',
                    'integrase', 'portal', 'tail', 'baseplate']):
                phage_genes.append((int(fields[3]), int(fields[4])))

    if not phage_genes:
        return []

    phage_genes.sort()

    # Cluster within 15kb
    clusters = []
    cur = [phage_genes[0]]
    for g in phage_genes[1:]:
        if g[0] - cur[-1][1] < 15000:
            cur.append(g)
        else:
            if len(cur) >= 2:
                clusters.append((cur[0][0], cur[-1][1]))
            cur = [g]
    if len(cur) >= 2:
        clusters.append((cur[0][0], cur[-1][1]))

    # Extend island boundaries by 2kb each side
    return [(max(0, s - 2000), e + 2000) for s, e in clusters]

def in_phage(start, end, islands):
    """Check if gene overlaps any phage island."""
    for (ps, pe) in islands:
        if start < pe and end > ps:
            return True
    return False

def match_preds(preds, refs):
    """Match predictions to reference, return per-prediction TP/FP status."""
    is_tp = [False] * len(preds)
    matched_ref = set()
    for pi, p in enumerate(preds):
        for ri, r in enumerate(refs):
            if ri in matched_ref: continue
            if p["strand"] != r["strand"]: continue
            ov_s = max(p["start"], r["start"])
            ov_e = min(p["end"], r["end"])
            if ov_s >= ov_e: continue
            ov = ov_e - ov_s
            if min(ov/max(p["end"]-p["start"],1), ov/max(r["end"]-r["start"],1)) >= 0.8:
                is_tp[pi] = True
                matched_ref.add(ri)
                break
    return is_tp, matched_ref

def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "/Users/akomissarov/Dropbox/workspace/new/renorm/data"
    pred_dir = sys.argv[2] if len(sys.argv) > 2 else "/tmp"

    genomes = [
        "ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
        "pseudomonas_aeruginosa", "saureus_nctc8325",
        "mycobacterium_tuberculosis", "synechocystis_pcc6803",
        "bacteroides_fragilis", "borrelia_burgdorferi"
    ]

    print(f"{'Genome':<25} {'Islands':>7} {'Size_kb':>7} | {'Phage':>6} {'Core':>6} | {'TP_ph':>5} {'FP_ph':>5} {'FN_ph':>5} | {'TP_co':>5} {'FP_co':>5} {'FN_co':>5} | {'F1_ph':>6} {'F1_co':>6}")
    print("-" * 120)

    for name in genomes:
        pred_path = os.path.join(pred_dir, f"{name}_test.gff")
        ref_path = os.path.join(data_dir, f"{name}.gff")
        if not os.path.exists(pred_path) or not os.path.exists(ref_path):
            continue

        preds = parse_prediction_gff(pred_path)
        refs = parse_reference_gff(ref_path)
        islands = find_phage_regions(ref_path)

        is_tp, matched_ref = match_preds(preds, refs)

        # Split by phage/core
        tp_ph, fp_ph, tp_co, fp_co = 0, 0, 0, 0
        for pi, p in enumerate(preds):
            phage = in_phage(p["start"], p["end"], islands)
            if is_tp[pi]:
                if phage: tp_ph += 1
                else: tp_co += 1
            else:
                if phage: fp_ph += 1
                else: fp_co += 1

        fn_ph, fn_co = 0, 0
        for ri, r in enumerate(refs):
            if ri not in matched_ref:
                if in_phage(r["start"], r["end"], islands):
                    fn_ph += 1
                else:
                    fn_co += 1

        # Count reference genes in phage regions
        ref_phage = sum(1 for r in refs if in_phage(r["start"], r["end"], islands))
        ref_core = len(refs) - ref_phage

        island_size = sum(e - s for s, e in islands) / 1000

        f1_ph = 2*tp_ph/(2*tp_ph+fp_ph+fn_ph) if (2*tp_ph+fp_ph+fn_ph) > 0 else 0
        f1_co = 2*tp_co/(2*tp_co+fp_co+fn_co) if (2*tp_co+fp_co+fn_co) > 0 else 0

        print(f"{name:<25} {len(islands):>7} {island_size:>7.0f} | {ref_phage:>6} {ref_core:>6} | {tp_ph:>5} {fp_ph:>5} {fn_ph:>5} | {tp_co:>5} {fp_co:>5} {fn_co:>5} | {f1_ph:>6.3f} {f1_co:>6.3f}")

if __name__ == "__main__":
    main()
