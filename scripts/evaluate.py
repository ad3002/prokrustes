#!/usr/bin/env python3
"""Evaluate GFF predictions against a reference GFF.

Usage:
    python evaluate.py prediction.gff reference.gff
    python evaluate.py prediction.gff reference.gff --tsv

Matching criterion: 80% reciprocal overlap on same strand.
Start accuracy: |predicted_start - known_start| <= 3 bp.
"""

import sys


def parse_prediction_gff(path):
    """Parse CDS entries from a prediction GFF (Prokrustes or Prodigal)."""
    genes = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'CDS':
                continue
            genes.append({
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
            })
    return genes


def parse_reference_gff(path):
    """Parse reference CDS from NCBI GFF3.

    We compare CDS to CDS — not gene to CDS — because:
    - Prokrustes predicts CDS (start codon to stop codon)
    - NCBI 'gene' entries may include UTRs and be longer
    - CDS entries are the actual coding sequences

    Deduplicates by (start, end, strand) to handle multiple
    protein products per locus. Skips pseudogenes.
    """
    cds_list = []
    seen = set()
    # Track parent gene pseudo status
    pseudo_genes = set()

    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            attrs = fields[8]

            # Collect pseudo gene IDs
            if fields[2] == 'gene':
                if 'pseudo=true' in attrs or 'pseudogene=' in attrs:
                    for part in attrs.split(';'):
                        if part.startswith('ID='):
                            pseudo_genes.add(part.split('=')[1])

            if fields[2] != 'CDS':
                continue

            # Skip CDS from pseudo genes
            is_pseudo = 'pseudo=true' in attrs or 'pseudogene=' in attrs
            if not is_pseudo:
                for part in attrs.split(';'):
                    if part.startswith('Parent='):
                        parent = part.split('=')[1]
                        # Parent can be gene-XXX or a list
                        for p in parent.split(','):
                            if p in pseudo_genes:
                                is_pseudo = True
                                break
            if is_pseudo:
                continue

            key = (int(fields[3]), int(fields[4]), fields[6])
            if key not in seen:
                seen.add(key)
                cds_list.append({
                    "start": key[0],
                    "end": key[1],
                    "strand": key[2],
                })

    return cds_list


def evaluate(predicted, known):
    """Compute TP, FP, FN, precision, recall, F1, start accuracy."""
    tp, fp = 0, 0
    correct_starts = 0
    matched = [False] * len(known)

    for pred in predicted:
        best_ov, best_idx = 0.0, None
        for i, gene in enumerate(known):
            if matched[i]:
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
                best_ov, best_idx = recip, i
        if best_ov >= 0.8 and best_idx is not None:
            tp += 1
            matched[best_idx] = True
            if abs(pred["start"] - known[best_idx]["start"]) <= 3:
                correct_starts += 1
        else:
            fp += 1

    fn = sum(1 for m in matched if not m)
    n_pred = len(predicted)
    n_known = len(known)
    sens = tp / n_known if n_known else 0
    prec = tp / n_pred if n_pred else 0
    f1 = 2 * sens * prec / (sens + prec) if (sens + prec) > 0 else 0
    start_acc = correct_starts / tp if tp > 0 else 0

    return {
        "predicted": n_pred, "known": n_known,
        "tp": tp, "fp": fp, "fn": fn,
        "sensitivity": sens, "precision": prec, "f1": f1,
        "start_accuracy": start_acc,
    }


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} prediction.gff reference.gff [--tsv]")
        sys.exit(1)

    pred_path = sys.argv[1]
    ref_path = sys.argv[2]
    tsv_mode = "--tsv" in sys.argv

    predicted = parse_prediction_gff(pred_path)
    known = parse_reference_gff(ref_path)

    m = evaluate(predicted, known)

    if tsv_mode:
        print(f"{m['known']}\t{m['predicted']}\t{m['tp']}\t{m['fp']}\t{m['fn']}\t"
              f"{m['precision']:.4f}\t{m['sensitivity']:.4f}\t{m['f1']:.4f}\t"
              f"{m['start_accuracy']:.4f}")
    else:
        print(f"  Reference:  {m['known']} CDS")
        print(f"  Predicted:  {m['predicted']} CDS")
        print(f"  TP: {m['tp']}  FP: {m['fp']}  FN: {m['fn']}")
        print(f"  Sensitivity: {m['sensitivity']*100:.1f}%")
        print(f"  Precision:   {m['precision']*100:.1f}%")
        print(f"  F1:          {m['f1']:.4f}")
        print(f"  Start acc:   {m['start_accuracy']*100:.1f}%")


if __name__ == "__main__":
    main()
