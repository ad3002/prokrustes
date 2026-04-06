#!/usr/bin/env python3
"""Compare detected prophage regions with known phage genes from NCBI."""
import sys, os

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
    # Cluster within 15kb
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

def parse_stderr_regions(stderr_path):
    """Parse prophage regions from Prokrustes stderr."""
    # This won't work — regions aren't printed individually.
    # Instead, parse from the GFF output looking for patterns.
    return []

def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    known = find_known_phage_clusters(f"{data_dir}/ecoli_k12.gff")
    print(f"Known phage islands: {len(known)}")
    total_known = 0
    for s, e in known:
        total_known += e - s
        print(f"  {s:>8}..{e:>8} ({(e-s)/1000:.0f}kb)")
    print(f"  Total: {total_known/1000:.0f}kb")

if __name__ == "__main__":
    main()
