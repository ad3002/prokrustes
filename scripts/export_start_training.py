#!/usr/bin/env python3
"""Export pairwise start training data for LightGBM ranker.
For each stop group with a wrong start, output feature vectors
for all candidate starts, with the correct one labeled."""

import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from evaluate import parse_reference_gff

def parse_dump_tsv(path):
    orfs = []
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 26: continue
            orfs.append({
                'start': int(parts[0]),
                'end': int(parts[1]),
                'strand': parts[2],
                'length': int(parts[3]),
                'start_codon': parts[4],
                'is_longest': 1 if parts[5] == 'true' else 0,
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
                'frame': int(parts[26]),
            })
    return orfs

def parse_pk_gff(path):
    genes = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS': continue
            genes.append({
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6],
            })
    return genes

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"
    out_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'start_training_v2.tsv')

    genomes = ["ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
               "saureus_nctc8325", "bacteroides_fragilis",
               "pseudomonas_aeruginosa", "mycobacterium_tuberculosis"]

    # Features for pairwise ranking
    feature_names = [
        'rbs', 'rbs_pwm', 'start_type', 'hex_avg', 'frame_bias', 'edge',
        'length', 'score', 'start_nn', 'start_ctx', 'upstream_at',
        'leaderless', 'gc3_bias', 'is_longest', 'hex_cov', 'viterbi_frac',
        # Relative features (computed per group)
        'delta_hex_to_longest', 'frac_of_longest', 'delta_rbs_to_best',
        'delta_score_to_best', 'length_rank', 'hex_rank',
    ]

    total_groups = 0
    total_correct = 0

    with open(out_path, 'w') as out:
        # Header: qid, label, features...
        out.write('qid\tlabel\t' + '\t'.join(feature_names) + '\n')

        qid = 0
        for name in genomes:
            dump_path = f"/tmp/{name}_all_orfs.tsv"
            pred_path = f"/tmp/{name}_test.gff"
            ref_path = f"{data_dir}/{name}.gff"
            if not os.path.exists(dump_path) or not os.path.exists(pred_path): continue

            preds = parse_pk_gff(pred_path)
            refs = parse_reference_gff(ref_path)
            orfs = parse_dump_tsv(dump_path)

            # Index ORFs by (strand, stop_group)
            groups = {}
            for o in orfs:
                key = (o['strand'], o['stop_group'])
                if key not in groups:
                    groups[key] = []
                groups[key].append(o)

            # For each predicted gene, find reference match (same end ±3bp)
            for pred in preds:
                best_ref = None
                for ref in refs:
                    if pred['strand'] != ref['strand']: continue
                    if abs(pred['end'] - ref['end']) <= 3:
                        best_ref = ref
                        break
                if best_ref is None: continue

                # Find the stop group
                group_key = None
                for key, members in groups.items():
                    if key[0] != pred['strand']: continue
                    if any(abs(m['end'] - pred['end']) <= 3 for m in members):
                        group_key = key
                        break
                if group_key is None: continue

                candidates = groups[group_key]
                if len(candidates) < 2: continue

                # Label: which candidate matches the reference start?
                has_correct = False
                for c in candidates:
                    if abs(c['start'] - best_ref['start']) <= 3:
                        has_correct = True
                        break
                if not has_correct: continue

                # Compute group-level stats for relative features
                max_len = max(c['length'] for c in candidates)
                max_rbs = max(c['rbs'] for c in candidates)
                max_hex = max(c['hex_avg'] for c in candidates)
                max_score = max(c['score'] for c in candidates)
                sorted_by_len = sorted(candidates, key=lambda c: -c['length'])
                sorted_by_hex = sorted(candidates, key=lambda c: -c['hex_avg'])

                qid += 1
                for c in candidates:
                    # Label: 1 if this is the correct start, 0 otherwise
                    label = 1 if abs(c['start'] - best_ref['start']) <= 3 else 0
                    if label == 1: total_correct += 1

                    # Start type as numeric
                    st = {'ATG': 1.0, 'GTG': 0.55, 'TTG': 0.35}.get(c['start_codon'], 0.2)

                    # Relative features
                    delta_hex = c['hex_avg'] - max_hex
                    frac_longest = c['length'] / max(max_len, 1)
                    delta_rbs = c['rbs'] - max_rbs
                    delta_score = c['score'] - max_score
                    len_rank = sorted_by_len.index(c) / max(len(candidates) - 1, 1)
                    hex_rank = sorted_by_hex.index(c) / max(len(candidates) - 1, 1)

                    features = [
                        c['rbs'], c['rbs_pwm'], st, c['hex_avg'], c['frame_bias'],
                        c['edge'], c['length'], c['score'], c['start_nn'], c['start_ctx'],
                        c['upstream_at'], c['leaderless'], c['gc3_bias'], c['is_longest'],
                        c['hex_cov'], c['viterbi_frac'],
                        delta_hex, frac_longest, delta_rbs, delta_score,
                        len_rank, hex_rank,
                    ]

                    out.write(f'{qid}\t{label}\t' + '\t'.join(f'{f:.6f}' for f in features) + '\n')

                total_groups += 1

    print(f"Exported: {total_groups} stop groups, {total_correct} correct labels")
    print(f"Output: {out_path}")
    print(f"Features: {len(feature_names)}")
    print(f"\nTrain with:")
    print(f"  python3 scripts/train_start_ranker.py")

if __name__ == "__main__":
    main()
