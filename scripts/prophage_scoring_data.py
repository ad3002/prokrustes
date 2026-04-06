#!/usr/bin/env python3
"""Collect data for prophage scoring law.
1. Longest low-hex run distribution in core vs phage
2. Gene density independence from hex
3. Island boundary behavior
4. Cryptic vs structured prophage characteristics"""

import sys, os, statistics
sys.path.insert(0, os.path.dirname(__file__))

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

def in_island(s, e, islands):
    for ps, pe in islands:
        if s < pe and e > ps: return True
    return False

def parse_pk(path):
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
    return sorted(genes, key=lambda g: g['start'])

def find_low_hex_runs(genes, threshold):
    """Find all runs of consecutive low-hex genes. Returns list of run lengths."""
    runs = []
    current_run = 0
    for g in genes:
        if g['hex'] < threshold:
            current_run += 1
        else:
            if current_run > 0:
                runs.append(current_run)
            current_run = 0
    if current_run > 0:
        runs.append(current_run)
    return runs

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    for name in ["ecoli_k12", "salmonella_lt2", "bacillus_subtilis", "saureus_nctc8325"]:
        pred_path = f"/tmp/{name}_test.gff"
        ref_path = f"{data_dir}/{name}.gff"
        if not os.path.exists(pred_path): continue

        islands = find_known_phage_islands(ref_path)
        genes = parse_pk(pred_path)
        if not islands or not genes: continue

        # Compute threshold (10th percentile of confident genes)
        confident_hex = sorted([g['hex'] for g in genes if g['length'] >= 600])
        if len(confident_hex) < 30: continue
        threshold = confident_hex[len(confident_hex) // 10]

        # Split genes into core and phage
        core_genes = [g for g in genes if not in_island(g['start'], g['end'], islands)]
        phage_genes = [g for g in genes if in_island(g['start'], g['end'], islands)]

        print(f"\n{'='*70}")
        print(f"  {name} (threshold={threshold:.3f})")
        print(f"{'='*70}")

        # === 1. Low-hex run distribution: core vs phage ===
        print(f"\n  1. LOW-HEX RUN LENGTHS (consecutive genes with hex < {threshold:.3f})")

        core_runs = find_low_hex_runs(core_genes, threshold)
        phage_runs = find_low_hex_runs(phage_genes, threshold)

        if core_runs:
            core_runs_sorted = sorted(core_runs, reverse=True)
            print(f"    CORE ({len(core_runs)} runs):")
            print(f"      max={core_runs_sorted[0]}, 99th pct={core_runs_sorted[max(0,len(core_runs_sorted)//100)]}")
            print(f"      95th pct={core_runs_sorted[len(core_runs_sorted)*5//100]}")
            print(f"      90th pct={core_runs_sorted[len(core_runs_sorted)*10//100]}")
            # Distribution of run lengths
            for k in [1, 2, 3, 4, 5, 6, 7, 8, 10, 15]:
                n = sum(1 for r in core_runs if r >= k)
                print(f"      runs >= {k:2d}: {n:4d}")

        if phage_runs:
            phage_runs_sorted = sorted(phage_runs, reverse=True)
            print(f"    PHAGE ({len(phage_runs)} runs):")
            print(f"      max={phage_runs_sorted[0]}")
            for k in [1, 2, 3, 4, 5, 6, 7, 8, 10, 15]:
                n = sum(1 for r in phage_runs if r >= k)
                print(f"      runs >= {k:2d}: {n:4d}")

        # === 2. Gene density independence from hex ===
        print(f"\n  2. GENE DENSITY vs HEX (are they independent signals?)")

        # Compute local density: for each gene, count neighbors within 5kb
        for label, glist in [("CORE", core_genes), ("PHAGE", phage_genes)]:
            densities = []
            for i, g in enumerate(glist):
                near = sum(1 for g2 in glist
                          if abs(g2['start'] - g['start']) < 5000 and g2 is not g)
                densities.append(near / 10.0)  # genes per kb (5kb window / 10 = per kb)

            low_hex = [d for g, d in zip(glist, densities) if g['hex'] < threshold]
            high_hex = [d for g, d in zip(glist, densities) if g['hex'] >= threshold]

            if low_hex and high_hex:
                print(f"    {label}: low-hex density={statistics.mean(low_hex):.2f}, "
                      f"high-hex density={statistics.mean(high_hex):.2f} genes/kb")

        # === 3. Island boundary behavior ===
        print(f"\n  3. ISLAND BOUNDARY BEHAVIOR")
        for s, e in islands:
            island_genes = [g for g in genes if g['start'] >= s and g['end'] <= e]
            if len(island_genes) < 5: continue

            # First 3, last 3, center genes
            first3 = island_genes[:3]
            last3 = island_genes[-3:]
            center = island_genes[len(island_genes)//3 : 2*len(island_genes)//3]

            first3_hex = [g['hex'] for g in first3]
            last3_hex = [g['hex'] for g in last3]
            center_hex = [g['hex'] for g in center]

            if center_hex:
                print(f"    Island {s}..{e} ({len(island_genes)} genes):")
                print(f"      boundary (first 3): hex={[f'{h:.2f}' for h in first3_hex]}")
                print(f"      center:             hex avg={statistics.mean(center_hex):.3f}")
                print(f"      boundary (last 3):  hex={[f'{h:.2f}' for h in last3_hex]}")

        # === 4. lowHexFraction per island ===
        print(f"\n  4. PER-ISLAND lowHexFraction")
        for s, e in islands:
            island_genes = [g for g in genes if g['start'] >= s and g['end'] <= e]
            if not island_genes: continue
            n_low = sum(1 for g in island_genes if g['hex'] < threshold)
            frac = n_low / len(island_genes)
            size_kb = (e - s) / 1000
            density = len(island_genes) / size_kb
            print(f"    {s:>8}..{e:>8} ({size_kb:.0f}kb): {n_low}/{len(island_genes)} low-hex "
                  f"({frac:.0%}), density={density:.1f} g/kb")

if __name__ == "__main__":
    main()
