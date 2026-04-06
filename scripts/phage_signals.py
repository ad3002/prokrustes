#!/usr/bin/env python3
"""Analyze what signals distinguish phage regions from core genome.
Uses our prediction features to understand what's different in phage islands."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from evaluate import parse_reference_gff

def find_phage_regions(gff_path):
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
                clusters.append((cur[0][0] - 2000, cur[-1][1] + 2000))
            cur = [g]
    if len(cur) >= 2:
        clusters.append((cur[0][0] - 2000, cur[-1][1] + 2000))
    return clusters

def in_phage(start, end, islands):
    for (ps, pe) in islands:
        if start < pe and end > ps: return True
    return False

def parse_pk_features(path):
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

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    for name in ["ecoli_k12", "salmonella_lt2", "bacillus_subtilis"]:
        pred_path = f"/tmp/{name}_test.gff"
        ref_path = f"{data_dir}/{name}.gff"
        if not os.path.exists(pred_path): continue

        islands = find_phage_regions(ref_path)
        preds = parse_pk_features(pred_path)

        phage_preds = [p for p in preds if in_phage(p['start'], p['end'], islands)]
        core_preds = [p for p in preds if not in_phage(p['start'], p['end'], islands)]

        if not phage_preds: continue

        import statistics
        print(f"=== {name} ({len(phage_preds)} phage, {len(core_preds)} core predictions) ===")

        for feature in ['length', 'hex', 'rbs', 'score']:
            ph_vals = [p[feature] for p in phage_preds]
            co_vals = [p[feature] for p in core_preds]
            ph_med = statistics.median(ph_vals)
            co_med = statistics.median(co_vals)
            print(f"  {feature:>8}: phage median={ph_med:.3f}, core median={co_med:.3f}")

        # Gene density: genes per kb
        for islands_data, label, genes in [
            (islands, "phage", phage_preds),
            ([], "core", core_preds)
        ]:
            if label == "phage":
                total_kb = sum(e-s for s,e in islands) / 1000
            else:
                total_kb = 4600  # rough
            density = len(genes) / max(total_kb, 0.1)
            avg_len = statistics.mean([g['length'] for g in genes]) if genes else 0
            print(f"  {label}: {len(genes)} genes in {total_kb:.0f}kb = {density:.1f} genes/kb, avg_len={avg_len:.0f}bp")
        print()

if __name__ == "__main__":
    main()
