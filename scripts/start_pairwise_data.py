#!/usr/bin/env python3
"""Collect pairwise start data: for each wrong-start gene, compare
features of our chosen start vs the correct start from reference.
This is the training data for the pairwise start ranker."""

import os, sys, statistics
sys.path.insert(0, os.path.dirname(__file__))
from evaluate import parse_prediction_gff, parse_reference_gff

def parse_pk_with_features(path):
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
                'start_type': attrs.get('start_type', '?'),
            })
    return genes

def parse_dump_tsv(path):
    """Parse --dump-all-orfs TSV."""
    orfs = []
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 25: continue
            orfs.append({
                'start': int(parts[0]),
                'end': int(parts[1]),
                'strand': parts[2],
                'length': int(parts[3]),
                'start_codon': parts[4],
                'is_longest': parts[5] == 'true',
                'hex_avg': float(parts[6]),
                'hex_total': float(parts[7]),
                'frame_bias': float(parts[8]),
                'hex_cov': float(parts[9]),
                'rbs': float(parts[10]),
                'rbs_pwm': float(parts[11]),
                'gc3': float(parts[12]),
                'upstream_at': float(parts[13]),
                'mono': float(parts[14]),
                'edge': float(parts[15]),
                'score': float(parts[16]),
                'leaderless': float(parts[19]),
                'start_ctx': float(parts[20]),
                'gc3_bias': float(parts[21]),
                'viterbi_frac': float(parts[23]),
                'start_nn': float(parts[24]),
                'stop_group': int(parts[25]),
            })
    return orfs

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    all_pairs = []  # (correct_features, wrong_features, delta_bp)

    for name in ["ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
                  "saureus_nctc8325", "bacteroides_fragilis",
                  "pseudomonas_aeruginosa", "mycobacterium_tuberculosis"]:

        dump_path = f"/tmp/{name}_all_orfs.tsv"
        pred_path = f"/tmp/{name}_test.gff"
        ref_path = f"{data_dir}/{name}.gff"

        if not os.path.exists(dump_path):
            print(f"SKIP {name}: need --dump-all-orfs")
            continue
        if not os.path.exists(pred_path): continue

        preds = parse_pk_with_features(pred_path)
        refs = parse_reference_gff(ref_path)
        orfs = parse_dump_tsv(dump_path)

        # Index ORFs by (end, strand) = stop group
        orf_by_end = {}
        for o in orfs:
            key = (o['end'], o['strand'])
            if key not in orf_by_end:
                orf_by_end[key] = []
            orf_by_end[key].append(o)

        n_pairs = 0

        for pred in preds:
            # Find matching reference (same end ±3bp)
            best_ref = None
            for ref in refs:
                if pred['strand'] != ref['strand']: continue
                if abs(pred['end'] - ref['end']) <= 3:
                    best_ref = ref
                    break
            if best_ref is None: continue

            start_diff = pred['start'] - best_ref['start']
            if abs(start_diff) <= 30: continue  # not a wrong start

            # Find both ORFs in dump (our chosen start and correct start)
            key = (pred['end'], pred['strand'])
            candidates = orf_by_end.get(key, [])
            # Also check ±3bp end
            for de in [-3, -2, -1, 0, 1, 2, 3]:
                k2 = (pred['end'] + de, pred['strand'])
                if k2 != key and k2 in orf_by_end:
                    candidates.extend(orf_by_end[k2])

            # Find our chosen ORF
            our_orf = None
            for c in candidates:
                if abs(c['start'] - pred['start']) <= 3:
                    our_orf = c
                    break

            # Find correct ORF
            correct_orf = None
            for c in candidates:
                if abs(c['start'] - best_ref['start']) <= 3:
                    correct_orf = c
                    break

            if our_orf is None or correct_orf is None: continue

            all_pairs.append({
                'genome': name,
                'our': our_orf,
                'correct': correct_orf,
                'delta_bp': start_diff,  # positive = we chose downstream (too short gene)
            })
            n_pairs += 1

        print(f"{name}: {n_pairs} wrong-start pairs found")

    print(f"\nTotal: {len(all_pairs)} pairs")

    if not all_pairs:
        print("No pairs found. Generate dumps first:")
        print("  target/release/prokrustes GENOME.fasta --dump-all-orfs > /tmp/NAME_all_orfs.tsv")
        return

    # Analyze feature deltas (correct - wrong)
    print(f"\n{'='*60}")
    print(f"  FEATURE DELTAS (correct_start - our_start)")
    print(f"  Positive = correct has MORE of this feature")
    print(f"{'='*60}")

    features = ['rbs', 'rbs_pwm', 'hex_avg', 'frame_bias', 'upstream_at',
                'leaderless', 'start_ctx', 'start_nn', 'edge', 'gc3_bias']

    for feat in features:
        deltas = [p['correct'][feat] - p['our'][feat] for p in all_pairs]
        if not deltas: continue
        med = statistics.median(deltas)
        pos = sum(1 for d in deltas if d > 0.01)
        neg = sum(1 for d in deltas if d < -0.01)
        print(f"  Δ{feat:>12}: median={med:+.4f}  correct>ours={pos} ({pos*100//len(deltas)}%)  ours>correct={neg} ({neg*100//len(deltas)}%)")

    # Length delta
    len_deltas = [p['correct']['length'] - p['our']['length'] for p in all_pairs]
    med = statistics.median(len_deltas)
    shorter = sum(1 for d in len_deltas if d < 0)
    print(f"\n  Δlength: median={med:+.0f}bp  correct shorter={shorter}/{len(all_pairs)} ({shorter*100//len(all_pairs)}%)")

    # Start codon type comparison
    our_atg = sum(1 for p in all_pairs if p['our']['start_codon'] == 'ATG')
    cor_atg = sum(1 for p in all_pairs if p['correct']['start_codon'] == 'ATG')
    print(f"\n  Start codon: our ATG={our_atg}/{len(all_pairs)}, correct ATG={cor_atg}/{len(all_pairs)}")

    # is_longest comparison
    our_long = sum(1 for p in all_pairs if p['our']['is_longest'])
    cor_long = sum(1 for p in all_pairs if p['correct']['is_longest'])
    print(f"  is_longest: ours={our_long}/{len(all_pairs)}, correct={cor_long}/{len(all_pairs)}")

    # too_long vs too_short breakdown
    too_long = [p for p in all_pairs if p['delta_bp'] < -30]  # our start upstream = too long
    too_short = [p for p in all_pairs if p['delta_bp'] > 30]   # our start downstream = too short
    print(f"\n  too_long (we chose upstream start): {len(too_long)}")
    print(f"  too_short (we chose downstream start): {len(too_short)}")

    # For too_long: what features does the correct (shorter) start have MORE of?
    if too_long:
        print(f"\n  TOO_LONG errors ({len(too_long)}):")
        print(f"    Our start is UPSTREAM → gene is longer than it should be")
        print(f"    Correct start has:")
        for feat in ['rbs', 'rbs_pwm', 'start_ctx', 'start_nn', 'upstream_at']:
            deltas = [p['correct'][feat] - p['our'][feat] for p in too_long]
            med = statistics.median(deltas)
            pos = sum(1 for d in deltas if d > 0.01)
            print(f"      Δ{feat}: median={med:+.4f}  correct>ours={pos}/{len(too_long)}")

    if too_short:
        print(f"\n  TOO_SHORT errors ({len(too_short)}):")
        print(f"    Our start is DOWNSTREAM → gene is shorter than it should be")
        for feat in ['rbs', 'rbs_pwm', 'start_ctx', 'start_nn', 'hex_avg']:
            deltas = [p['correct'][feat] - p['our'][feat] for p in too_short]
            med = statistics.median(deltas)
            pos = sum(1 for d in deltas if d > 0.01)
            print(f"      Δ{feat}: median={med:+.4f}  correct>ours={pos}/{len(too_short)}")

if __name__ == "__main__":
    main()
