#!/usr/bin/env python3
"""Validate prophage detection: compare detected islands with known phage clusters.
Measures island-level precision, recall, and boundary accuracy."""

import sys, os, re

def find_known_phage_clusters(gff_path):
    """Extract known phage gene clusters from NCBI GFF."""
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
                clusters.append((cur[0][0], cur[-1][1], len(cur)))
            cur = [g]
    if len(cur) >= 2:
        clusters.append((cur[0][0], cur[-1][1], len(cur)))
    return clusters

def parse_detected_regions(stderr_text):
    """Parse detected regions from Prokrustes stderr."""
    regions = []
    for line in stderr_text.split('\n'):
        m = re.match(r'\s*Region:\s+(\d+)\.\.(\d+)\s+\([\d.]+kb\)\s+(\w+)', line)
        if m:
            regions.append((int(m.group(1)), int(m.group(2)), m.group(3)))
    return regions

def overlap_fraction(s1, e1, s2, e2):
    """Fraction of (s1,e1) overlapped by (s2,e2)."""
    ov = max(0, min(e1, e2) - max(s1, s2))
    return ov / max(e1 - s1, 1)

def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "/Users/akomissarov/Dropbox/workspace/new/renorm/data"
    stderr_dir = sys.argv[2] if len(sys.argv) > 2 else "/tmp"

    genomes = [
        "ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
        "saureus_nctc8325", "bacteroides_fragilis", "borrelia_burgdorferi",
    ]

    for name in genomes:
        ref_path = os.path.join(data_dir, f"{name}.gff")
        stderr_path = os.path.join(stderr_dir, f"{name}_stderr.txt")

        if not os.path.exists(ref_path):
            continue

        known = find_known_phage_clusters(ref_path)
        if not known:
            print(f"{name}: no known phage clusters")
            continue

        # Run Prokrustes and capture stderr
        import subprocess
        fasta = os.path.join(data_dir, f"{name}.fasta")
        ncrna = os.path.join(data_dir, f"{name}_ncrna.gff")
        ncrna_flag = f"--ncrna {ncrna}" if os.path.exists(ncrna) else ""

        binary = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                              "target/release/prokrustes")
        cmd = f"{binary} {fasta} {ncrna_flag}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)

        detected = parse_detected_regions(result.stderr)

        print(f"\n{'='*60}")
        print(f"  {name}")
        print(f"{'='*60}")
        print(f"\nKnown phage islands ({len(known)}):")
        total_known_kb = 0
        for s, e, n_genes in known:
            total_known_kb += (e - s) / 1000
            # Check if detected
            best_ov = 0
            best_det = None
            for ds, de, dt in detected:
                ov = overlap_fraction(s, e, ds, de)
                if ov > best_ov:
                    best_ov = ov
                    best_det = (ds, de, dt)
            status = f"DETECTED ({best_ov*100:.0f}% overlap, {best_det[2]})" if best_ov > 0.3 else "MISSED"
            print(f"  {s:>8}..{e:>8} ({(e-s)/1000:.0f}kb, {n_genes} genes) → {status}")
        print(f"  Total known: {total_known_kb:.0f}kb")

        print(f"\nDetected regions ({len(detected)}):")
        total_det_kb = 0
        for ds, de, dt in detected:
            total_det_kb += (de - ds) / 1000
            # Check if matches any known
            best_ov = 0
            best_known = None
            for s, e, ng in known:
                ov = overlap_fraction(ds, de, s, e)
                if ov > best_ov:
                    best_ov = ov
                    best_known = (s, e)
            if best_ov > 0.3:
                status = f"TRUE POSITIVE (matches {best_known[0]}..{best_known[1]})"
            else:
                status = "FALSE POSITIVE"
            print(f"  {ds:>8}..{de:>8} ({(de-ds)/1000:.0f}kb) {dt:>10} → {status}")
        print(f"  Total detected: {total_det_kb:.0f}kb")

        # Summary metrics
        tp = sum(1 for s, e, _ in known
                 if any(overlap_fraction(s, e, ds, de) > 0.3 for ds, de, _ in detected))
        fn = len(known) - tp
        fp = sum(1 for ds, de, _ in detected
                 if not any(overlap_fraction(ds, de, s, e) > 0.3 for s, e, _ in known))

        prec = tp / max(tp + fp, 1)
        rec = tp / max(tp + fn, 1)
        f1 = 2 * prec * rec / max(prec + rec, 0.001)

        print(f"\n  Island-level: TP={tp} FP={fp} FN={fn}")
        print(f"  Precision={prec:.2f} Recall={rec:.2f} F1={f1:.2f}")

if __name__ == "__main__":
    main()
