#!/usr/bin/env python3
"""Learn gene properties from 10 reference genomes.
Extract distributions of: intergenic distances, gene lengths, overlap patterns,
start codon types, strand patterns. These should inform our pipeline parameters."""

import sys
import os
import statistics

def parse_cds(gff_path):
    """Parse CDS from NCBI GFF."""
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9: continue
            if fields[2] != 'CDS': continue
            attrs = fields[8]
            if 'pseudo=true' in attrs or 'pseudogene=' in attrs: continue
            genes.append({
                "seqid": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "length": int(fields[4]) - int(fields[3]) + 1,
            })
    # Deduplicate
    seen = set()
    unique = []
    for g in genes:
        key = (g["start"], g["end"], g["strand"])
        if key not in seen:
            seen.add(key)
            unique.append(g)
    return unique

def analyze_genome(name, gff_path):
    """Analyze one genome's gene properties."""
    genes = parse_cds(gff_path)
    if not genes:
        return None

    # Sort by start position
    genes.sort(key=lambda g: (g["seqid"], g["start"]))

    lengths = [g["length"] for g in genes]

    # Intergenic distances (same strand)
    same_strand_gaps = []
    opp_strand_gaps = []
    same_strand_overlaps = []
    opp_strand_overlaps = []

    for i in range(len(genes) - 1):
        a, b = genes[i], genes[i+1]
        if a["seqid"] != b["seqid"]: continue

        gap = b["start"] - a["end"]  # negative = overlap
        if a["strand"] == b["strand"]:
            if gap >= 0:
                same_strand_gaps.append(gap)
            else:
                same_strand_overlaps.append(-gap)
        else:
            if gap >= 0:
                opp_strand_gaps.append(gap)
            else:
                opp_strand_overlaps.append(-gap)

    # Length distribution
    short_genes = sum(1 for l in lengths if l < 300)
    medium_genes = sum(1 for l in lengths if 300 <= l < 600)
    long_genes = sum(1 for l in lengths if l >= 600)

    # Gene density
    total_coding = sum(lengths)
    # Rough genome size from last gene end
    genome_size = max(g["end"] for g in genes)
    density = total_coding / genome_size if genome_size > 0 else 0

    return {
        "name": name,
        "n_genes": len(genes),
        "genome_size": genome_size,
        "density": density,
        "lengths": lengths,
        "short_pct": short_genes / len(genes) * 100,
        "medium_pct": medium_genes / len(genes) * 100,
        "long_pct": long_genes / len(genes) * 100,
        "median_length": statistics.median(lengths),
        "same_strand_gaps": same_strand_gaps,
        "opp_strand_gaps": opp_strand_gaps,
        "same_strand_overlaps": same_strand_overlaps,
        "opp_strand_overlaps": opp_strand_overlaps,
    }

def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    genomes = [
        "ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
        "pseudomonas_aeruginosa", "saureus_nctc8325",
        "mycobacterium_tuberculosis", "mycoplasma_genitalium",
        "synechocystis_pcc6803", "bacteroides_fragilis",
        "borrelia_burgdorferi"
    ]

    all_gaps_ss = []
    all_gaps_os = []
    all_ov_ss = []
    all_ov_os = []
    all_lengths = []

    print("=" * 80)
    print("GENE PROPERTIES ACROSS 10 REFERENCE GENOMES")
    print("=" * 80)
    print()

    for name in genomes:
        gff = os.path.join(data_dir, f"{name}.gff")
        if not os.path.exists(gff):
            print(f"{name}: GFF not found")
            continue

        r = analyze_genome(name, gff)
        if not r: continue

        print(f"--- {name} ---")
        print(f"  Genes: {r['n_genes']}, Density: {r['density']:.1%}")
        print(f"  Length distribution: short(<300)={r['short_pct']:.0f}% medium(300-600)={r['medium_pct']:.0f}% long(>600)={r['long_pct']:.0f}%")
        print(f"  Median gene length: {r['median_length']:.0f}bp")

        if r['same_strand_gaps']:
            print(f"  Same-strand gaps: n={len(r['same_strand_gaps'])}, median={statistics.median(r['same_strand_gaps']):.0f}bp")
            # Count ATGA-like overlaps (-1 to -4)
            coupling = sum(1 for ov in r['same_strand_overlaps'] if 1 <= ov <= 4)
            print(f"  Same-strand overlaps: n={len(r['same_strand_overlaps'])}, coupling(1-4bp)={coupling}")

        if r['opp_strand_overlaps']:
            print(f"  Opp-strand overlaps: n={len(r['opp_strand_overlaps'])}, median={statistics.median(r['opp_strand_overlaps']):.0f}bp, max={max(r['opp_strand_overlaps'])}bp")

        all_gaps_ss.extend(r['same_strand_gaps'])
        all_gaps_os.extend(r['opp_strand_gaps'])
        all_ov_ss.extend(r['same_strand_overlaps'])
        all_ov_os.extend(r['opp_strand_overlaps'])
        all_lengths.extend(r['lengths'])
        print()

    print("=" * 80)
    print("AGGREGATE STATISTICS")
    print("=" * 80)
    print()

    if all_gaps_ss:
        print(f"Same-strand intergenic gaps (n={len(all_gaps_ss)}):")
        pcts = [0, 10, 25, 50, 75, 90, 95, 99, 100]
        vals = [sorted(all_gaps_ss)[int(len(all_gaps_ss)*p/100)] if p < 100 else max(all_gaps_ss) for p in pcts]
        print(f"  Percentiles: {dict(zip(pcts, vals))}")
        print(f"  Genes with gap=0 (adjacent): {sum(1 for g in all_gaps_ss if g == 0)}")
        print(f"  Genes with gap<50: {sum(1 for g in all_gaps_ss if g < 50)} ({sum(1 for g in all_gaps_ss if g < 50)*100//len(all_gaps_ss)}%)")
        print()

    if all_ov_ss:
        print(f"Same-strand overlaps (n={len(all_ov_ss)}):")
        print(f"  1-4bp (coupling): {sum(1 for o in all_ov_ss if 1<=o<=4)} ({sum(1 for o in all_ov_ss if 1<=o<=4)*100//len(all_ov_ss)}%)")
        print(f"  5-30bp: {sum(1 for o in all_ov_ss if 5<=o<=30)}")
        print(f"  >30bp: {sum(1 for o in all_ov_ss if o>30)}")
        print()

    if all_ov_os:
        print(f"Opposite-strand overlaps (n={len(all_ov_os)}):")
        pcts = [50, 75, 90, 95, 99, 100]
        vals = [sorted(all_ov_os)[int(len(all_ov_os)*p/100)] if p < 100 else max(all_ov_os) for p in pcts]
        print(f"  Percentiles: {dict(zip(pcts, vals))}")
        print(f"  >200bp: {sum(1 for o in all_ov_os if o>200)}")
        print()

    if all_lengths:
        print(f"Gene lengths (n={len(all_lengths)}):")
        pcts = [1, 5, 10, 25, 50, 75, 90, 95, 99]
        vals = [sorted(all_lengths)[int(len(all_lengths)*p/100)] for p in pcts]
        print(f"  Percentiles: {dict(zip(pcts, vals))}")
        short = sum(1 for l in all_lengths if l < 150)
        print(f"  <150bp: {short} ({short*100//len(all_lengths)}%)")
        tiny = sum(1 for l in all_lengths if l < 90)
        print(f"  <90bp: {tiny} ({tiny*100//len(all_lengths)}%)")

if __name__ == "__main__":
    main()
