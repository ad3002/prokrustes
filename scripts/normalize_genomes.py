#!/usr/bin/env python3
"""Normalize genomes: rotate to start at dnaA, flip if dnaA on (-) strand.
Output: normalized FASTA + updated GFF with recalculated coordinates."""

import os, sys

def read_fasta(path):
    records = []
    name = ""
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name: records.append((name, ''.join(parts)))
                name = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip().upper())
    if name: records.append((name, ''.join(parts)))
    return records

def revcomp(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp.get(c, 'N') for c in reversed(seq))

def find_dnaa(gff_path):
    """Find dnaA gene position and strand."""
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9: continue
            if fields[2] != 'gene': continue
            if 'gene=dnaA' not in fields[8]: continue
            return {
                'seqid': fields[0],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6],
            }
    return None

def normalize_genome(seq, dnaa_start, dnaa_strand, genome_len):
    """Rotate genome so dnaA starts at position 1, flip if on (-) strand."""
    # Step 1: rotate so dnaA start is at position 0
    rotate_by = dnaa_start - 1  # 0-based
    rotated = seq[rotate_by:] + seq[:rotate_by]

    # Step 2: if dnaA was on (-) strand, reverse complement
    flipped = False
    if dnaa_strand == '-':
        # After rotation, dnaA end is near position 0
        # After revcomp, dnaA will be at the start on (+) strand
        rotated = revcomp(rotated)
        flipped = True

    return rotated, rotate_by, flipped

def transform_coord(start, end, strand, rotate_by, flipped, genome_len):
    """Transform gene coordinates after rotation + optional flip."""
    # Step 1: rotate
    new_start = (start - 1 - rotate_by) % genome_len + 1
    new_end = (end - 1 - rotate_by) % genome_len + 1

    # Handle wrap-around (gene spans origin after rotation)
    if new_start > new_end:
        # Gene wraps — skip for now (rare)
        return None, None, None

    new_strand = strand

    # Step 2: flip
    if flipped:
        old_start = new_start
        old_end = new_end
        new_start = genome_len - old_end + 1
        new_end = genome_len - old_start + 1
        new_strand = '-' if strand == '+' else '+'

    return new_start, new_end, new_strand

def normalize_gff(gff_path, out_path, target_seqid, rotate_by, flipped, genome_len):
    """Rewrite GFF with transformed coordinates."""
    n_transformed = 0
    n_skipped = 0

    with open(gff_path) as f, open(out_path, 'w') as out:
        for line in f:
            if line.startswith('#'):
                out.write(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                out.write(line)
                continue

            if fields[0] != target_seqid:
                # Different contig — skip or pass through
                out.write(line)
                continue

            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            new_start, new_end, new_strand = transform_coord(
                start, end, strand, rotate_by, flipped, genome_len)

            if new_start is None:
                n_skipped += 1
                continue

            fields[3] = str(new_start)
            fields[4] = str(new_end)
            fields[6] = new_strand
            out.write('\t'.join(fields) + '\n')
            n_transformed += 1

    return n_transformed, n_skipped

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"
    out_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'normalized')
    os.makedirs(out_dir, exist_ok=True)

    genomes = [
        "ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
        "pseudomonas_aeruginosa", "saureus_nctc8325",
        "mycobacterium_tuberculosis", "mycoplasma_genitalium",
        "synechocystis_pcc6803", "bacteroides_fragilis",
        "borrelia_burgdorferi"
    ]

    print(f"{'Genome':<25} {'dnaA pos':>10} {'strand':>6} {'rotate':>8} {'flip':>5} {'new dnaA':>10}")
    print("-" * 70)

    for name in genomes:
        fasta_path = f"{data_dir}/{name}.fasta"
        gff_path = f"{data_dir}/{name}.gff"
        if not os.path.exists(fasta_path) or not os.path.exists(gff_path):
            print(f"{name}: MISSING")
            continue

        dnaa = find_dnaa(gff_path)
        if dnaa is None:
            print(f"{name}: dnaA NOT FOUND — skipping")
            continue

        # Read main chromosome (first/largest record)
        records = read_fasta(fasta_path)
        main_record = max(records, key=lambda r: len(r[1]))
        main_name, main_seq = main_record
        main_seqid = main_name.split()[0]
        genome_len = len(main_seq)

        # Check if dnaA is on the main chromosome
        # GFF seqid format varies — extract the accession
        gff_seqid = dnaa['seqid']

        # Normalize
        normalized_seq, rotate_by, flipped = normalize_genome(
            main_seq, dnaa['start'], dnaa['strand'], genome_len)

        # Verify dnaA is now at start
        new_start, new_end, new_strand = transform_coord(
            dnaa['start'], dnaa['end'], dnaa['strand'],
            rotate_by, flipped, genome_len)

        flip_str = "YES" if flipped else "no"
        print(f"{name:<25} {dnaa['start']:>10} {dnaa['strand']:>6} {rotate_by:>8} {flip_str:>5} {new_start}..{new_end} ({new_strand})")

        # Write normalized FASTA
        out_fasta = os.path.join(out_dir, f"{name}.fasta")
        with open(out_fasta, 'w') as f:
            f.write(f">{main_seqid} normalized dnaA=1\n")
            for i in range(0, len(normalized_seq), 80):
                f.write(normalized_seq[i:i+80] + '\n')

        # Write normalized GFF
        out_gff = os.path.join(out_dir, f"{name}.gff")
        n_ok, n_skip = normalize_gff(gff_path, out_gff, gff_seqid,
                                      rotate_by, flipped, genome_len)
        print(f"  → {n_ok} features transformed, {n_skip} skipped (wrap-around)")

    # Verify: check strand distribution after normalization
    print(f"\n{'='*60}")
    print(f"  STRAND CHECK AFTER NORMALIZATION")
    print(f"{'='*60}")

    def get_gene_strands(gff_path):
        genes = {}
        with open(gff_path) as f:
            for line in f:
                if line.startswith('#'): continue
                fields = line.strip().split('\t')
                if len(fields) < 9 or fields[2] != 'gene': continue
                if 'gene_biotype=protein_coding' not in fields[8]: continue
                for part in fields[8].split(';'):
                    if part.startswith('gene='):
                        genes[part.split('=')[1]] = fields[6]
                        break
        return genes

    ecoli_file = os.path.join(out_dir, "ecoli_k12.gff")
    if os.path.exists(ecoli_file):
        ecoli_strands = get_gene_strands(ecoli_file)

        for name in genomes[1:]:
            other_file = os.path.join(out_dir, f"{name}.gff")
            if not os.path.exists(other_file): continue
            other_strands = get_gene_strands(other_file)

            same = sum(1 for g in ecoli_strands if g in other_strands and ecoli_strands[g] == other_strands[g])
            total = sum(1 for g in ecoli_strands if g in other_strands)
            pct = same * 100 // max(total, 1)
            print(f"  ecoli ↔ {name:<25}: {same}/{total} same strand ({pct}%)")

if __name__ == "__main__":
    main()
