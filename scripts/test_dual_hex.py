#!/usr/bin/env python3
"""Test dual hex model idea: would a local phage hex model
give better scores to phage genes than the host model?

For each known phage island:
1. Collect ORFs inside (from --dump-all-orfs)
2. Show their host hex scores
3. Estimate what a local model would give

Key question: are phage genes low-hex because wrong model,
or because genuinely weak coding signal?"""

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

def parse_dump_tsv(path):
    """Parse --dump-all-orfs TSV output."""
    orfs = []
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 17: continue
            orfs.append({
                'start': int(parts[0]),
                'end': int(parts[1]),
                'strand': parts[2],
                'length': int(parts[3]),
                'hex_avg': float(parts[6]),
                'frame_bias': float(parts[8]),
                'rbs': float(parts[10]),
                'score': float(parts[16]),
                'is_longest': parts[5] == 'true',
            })
    return orfs

def in_island(s, e, islands):
    for ps, pe in islands:
        if s >= ps and e <= pe: return True
    return False

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    for name in ["ecoli_k12", "salmonella_lt2"]:
        dump_path = f"/tmp/{name}_all_orfs.tsv"
        ref_path = f"{data_dir}/{name}.gff"

        if not os.path.exists(dump_path):
            print(f"Run: target/release/prokrustes {data_dir}/{name}.fasta --dump-all-orfs > {dump_path}")
            continue

        islands = find_known_phage_islands(ref_path)
        orfs = parse_dump_tsv(dump_path)
        if not islands or not orfs: continue

        print(f"\n{'='*70}")
        print(f"  {name}: {len(islands)} islands, {len(orfs)} total ORFs")
        print(f"{'='*70}")

        # Split ORFs by location
        phage_orfs = [o for o in orfs if in_island(o['start'], o['end'], islands)]
        core_orfs = [o for o in orfs if not in_island(o['start'], o['end'], islands)]

        # Compare hex distributions
        phage_hex = [o['hex_avg'] for o in phage_orfs if o['length'] >= 300]
        core_hex = [o['hex_avg'] for o in core_orfs if o['length'] >= 300]

        if phage_hex and core_hex:
            print(f"\n  Hex scores (ORFs ≥300bp):")
            print(f"    Core ({len(core_hex)}):  median={statistics.median(core_hex):.3f}")
            print(f"    Phage ({len(phage_hex)}): median={statistics.median(phage_hex):.3f}")

        # Key question: within a phage island, do ORFs score each other better?
        # If we train hex on phage genes and score phage genes with it,
        # do they look "coding" or still "noncoding"?

        # Proxy: frame_bias. If phage genes have positive frame bias,
        # they ARE coding (just different codon usage from host).
        # If frame_bias ≈ 0, they're genuinely ambiguous.
        phage_fb = [o['frame_bias'] for o in phage_orfs if o['length'] >= 300]
        core_fb = [o['frame_bias'] for o in core_orfs if o['length'] >= 300]

        if phage_fb and core_fb:
            print(f"\n  Frame bias (codon asymmetry — positive = coding):")
            print(f"    Core:  median={statistics.median(core_fb):.3f}")
            print(f"    Phage: median={statistics.median(phage_fb):.3f}")
            phage_pos_fb = sum(1 for fb in phage_fb if fb > 0.1)
            print(f"    Phage with positive fb: {phage_pos_fb}/{len(phage_fb)} ({phage_pos_fb*100//len(phage_fb)}%)")

        # Per-island analysis
        print(f"\n  Per-island hex and frame_bias:")
        for s, e in islands:
            island_orfs = [o for o in orfs if o['start'] >= s and o['end'] <= e and o['length'] >= 200]
            if len(island_orfs) < 3: continue

            hexes = [o['hex_avg'] for o in island_orfs]
            fbs = [o['frame_bias'] for o in island_orfs]
            n_pos_fb = sum(1 for fb in fbs if fb > 0.1)

            # If a local hex model were trained on THIS island's ORFs,
            # they would all score ~0 (by definition — they ARE the model).
            # The question is: would they score HIGHER than with host model?
            # Answer: if frame_bias is positive, yes — the codon usage is
            # different from host but internally consistent.

            print(f"    {s:>8}..{e:>8} ({(e-s)//1000}kb, {len(island_orfs)} ORFs≥200): "
                  f"hex={statistics.median(hexes):.3f}, fb={statistics.median(fbs):.3f}, "
                  f"pos_fb={n_pos_fb}/{len(island_orfs)}")

        # Estimate improvement from dual model
        # Phage genes with low hex but positive frame_bias = would benefit from local model
        beneficiaries = [o for o in phage_orfs if o['length'] >= 200
                        and o['hex_avg'] < statistics.median(core_hex) * 0.5  # low host hex
                        and o['frame_bias'] > 0.05]  # but coding signal present

        total_phage_200 = len([o for o in phage_orfs if o['length'] >= 200])
        print(f"\n  Potential beneficiaries (low host-hex, positive frame_bias):")
        print(f"    {len(beneficiaries)} / {total_phage_200} phage ORFs ≥200bp")
        print(f"    These would score higher with a local phage hex model")

if __name__ == "__main__":
    main()
