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

    NCBI GFF has multiple CDS rows per gene (one per protein/exon).
    We deduplicate by taking unique (start, end, strand) tuples.
    We also skip pseudogenes and non-protein-coding features.
    """
    genes = []
    seen = set()
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            # Use 'gene' type for NCBI reference — one per locus
            # Fall back to 'CDS' if no 'gene' entries found
            if fields[2] not in ('gene', 'CDS'):
                continue
            attrs = fields[8]
            # Skip pseudogenes
            if 'pseudo=true' in attrs or 'pseudogene=' in attrs:
                continue
            # Skip non-protein-coding genes (rRNA, tRNA, ncRNA)
            if 'gene_biotype=' in attrs and 'protein_coding' not in attrs:
                # Only skip if biotype is explicitly non-coding
                biotype = ''
                for part in attrs.split(';'):
                    if part.startswith('gene_biotype='):
                        biotype = part.split('=')[1]
                if biotype and biotype != 'protein_coding':
                    continue

            key = (int(fields[3]), int(fields[4]), fields[6])
            if key not in seen:
                seen.add(key)
                genes.append({
                    "start": key[0],
                    "end": key[1],
                    "strand": key[2],
                })
    return genes


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
