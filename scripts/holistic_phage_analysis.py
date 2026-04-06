#!/usr/bin/env python3
"""Holistic analysis: use ALL available signals to understand phage regions.
Gene-level hex scores, tRNA positions, IS elements, gene density, sizes."""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from evaluate import parse_reference_gff

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

def find_trna_positions(gff_path):
    trnas = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9: continue
            if fields[2] == 'tRNA':
                trnas.append((int(fields[3]), int(fields[4])))
    return sorted(trnas)

def parse_pk_genes(path):
    genes = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS': continue
            attrs = {}
            for part in fields[8].split(';'):
                if '=' in part: k, v = part.split('=', 1); attrs[k] = v
            genes.append({
                'start': int(fields[3]), 'end': int(fields[4]),
                'strand': fields[6],
                'length': int(fields[4]) - int(fields[3]) + 1,
                'hex': float(attrs.get('hex_score', '0')),
                'rbs': float(attrs.get('rbs_score', '0')),
                'score': float(fields[5]),
            })
    return genes

def in_island(pos_s, pos_e, islands):
    for s, e in islands:
        if pos_s < e and pos_e > s: return True
    return False

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    for name in ["ecoli_k12", "salmonella_lt2", "bacillus_subtilis", "saureus_nctc8325"]:
        ref_path = f"{data_dir}/{name}.gff"
        pred_path = f"/tmp/{name}_test.gff"
        if not os.path.exists(pred_path): continue

        islands = find_known_phage_islands(ref_path)
        trnas = find_trna_positions(ref_path)
        preds = parse_pk_genes(pred_path)

        if not islands: continue

        print(f"\n{'='*70}")
        print(f"  {name}")
        print(f"{'='*70}")

        # 1. Gene-level hex scores: phage vs core
        phage_hex = [g['hex'] for g in preds if in_island(g['start'], g['end'], islands)]
        core_hex = [g['hex'] for g in preds if not in_island(g['start'], g['end'], islands)]

        import statistics
        if phage_hex and core_hex:
            print(f"\n  Gene hex_avg:")
            print(f"    Core:  median={statistics.median(core_hex):.3f}, mean={statistics.mean(core_hex):.3f}")
            print(f"    Phage: median={statistics.median(phage_hex):.3f}, mean={statistics.mean(phage_hex):.3f}")
            # What fraction of phage genes have hex < core 10th percentile?
            core_p10 = sorted(core_hex)[len(core_hex)//10]
            low_hex_phage = sum(1 for h in phage_hex if h < core_p10)
            print(f"    Core 10th pct: {core_p10:.3f}")
            print(f"    Phage genes below core-10th: {low_hex_phage}/{len(phage_hex)} ({low_hex_phage*100//len(phage_hex)}%)")

        # 2. tRNA proximity to phage islands
        print(f"\n  tRNA near phage islands:")
        for s, e in islands:
            nearby_trna = [(ts, te) for ts, te in trnas
                          if abs(ts - s) < 5000 or abs(te - e) < 5000
                          or abs(ts - e) < 5000 or abs(te - s) < 5000]
            size_kb = (e - s) / 1000
            has_trna = "YES" if nearby_trna else "no"
            # Count our predicted genes in this island
            n_pred = sum(1 for g in preds if in_island(g['start'], g['end'], [(s, e)]))
            # Average hex of our predictions here
            island_hex = [g['hex'] for g in preds if in_island(g['start'], g['end'], [(s, e)])]
            avg_hex = statistics.mean(island_hex) if island_hex else 0
            print(f"    {s:>8}..{e:>8} ({size_kb:.0f}kb) tRNA_adjacent={has_trna} "
                  f"pred_genes={n_pred} avg_hex={avg_hex:.3f}")

        # 3. Island size distribution
        sizes = [(e-s)/1000 for s, e in islands]
        print(f"\n  Island sizes: {[f'{s:.0f}kb' for s in sizes]}")
        print(f"    Range: {min(sizes):.0f}-{max(sizes):.0f}kb, median={statistics.median(sizes):.0f}kb")

        # 4. Gene density in islands vs core
        total_island_bp = sum(e-s for s, e in islands)
        total_core_bp = max(g['end'] for g in preds) - total_island_bp
        n_island_genes = len(phage_hex)
        n_core_genes = len(core_hex)
        island_density = n_island_genes / (total_island_bp/1000) if total_island_bp > 0 else 0
        core_density = n_core_genes / (total_core_bp/1000) if total_core_bp > 0 else 0
        print(f"\n  Gene density:")
        print(f"    Core:  {core_density:.2f} genes/kb")
        print(f"    Phage: {island_density:.2f} genes/kb")

        # 5. Consecutive low-hex gene runs (potential outlier clusters)
        print(f"\n  Low-hex gene clusters (≥3 consecutive genes with hex < core 10th pct):")
        if core_hex:
            core_p10 = sorted(core_hex)[len(core_hex)//10]
            sorted_preds = sorted(preds, key=lambda g: g['start'])
            run_start = None
            run_count = 0
            for g in sorted_preds:
                if g['hex'] < core_p10:
                    if run_start is None:
                        run_start = g['start']
                    run_count += 1
                else:
                    if run_count >= 3:
                        run_end = sorted_preds[sorted_preds.index(g) - 1]['end'] if g != sorted_preds[0] else g['end']
                        is_real = in_island(run_start, run_end, islands)
                        print(f"    {run_start:>8}..{run_end:>8} ({run_count} genes) "
                              f"{'REAL PHAGE' if is_real else 'false alarm'}")
                    run_start = None
                    run_count = 0
            if run_count >= 3:
                print(f"    {run_start:>8}..{sorted_preds[-1]['end']:>8} ({run_count} genes)")

if __name__ == "__main__":
    main()
