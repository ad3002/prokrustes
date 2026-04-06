#!/usr/bin/env python3
"""Build integrase detector from reference annotations.
Extract DNA sequences of all integrases, find characteristic k-mer patterns."""

import sys, os

def extract_marker_genes(gff_path, fasta_path, keywords):
    """Extract DNA sequences of genes matching keywords."""
    # Read FASTA
    seqs = {}
    name = ""
    parts = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                if name: seqs[name] = ''.join(parts)
                name = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip().upper())
    if name: seqs[name] = ''.join(parts)

    # Extract marker gene sequences
    genes = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS': continue
            attrs = fields[8].lower()
            if not any(kw in attrs for kw in keywords): continue
            if 'pseudo=true' in attrs: continue

            seqid = fields[0]
            start = int(fields[3]) - 1  # 0-based
            end = int(fields[4])
            strand = fields[6]

            if seqid not in seqs: continue
            seq = seqs[seqid][start:end]
            if strand == '-':
                comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                seq = ''.join(comp.get(c, 'N') for c in reversed(seq))

            # Get product name
            product = ""
            for part in fields[8].split(';'):
                if part.startswith('product='):
                    product = part.split('=', 1)[1]

            genes.append({
                'seqid': seqid, 'start': start, 'end': end,
                'strand': strand, 'seq': seq, 'length': len(seq),
                'product': product,
            })
    return genes

def kmer_profile(seq, k=6):
    """Compute k-mer frequency profile."""
    counts = {}
    total = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if 'N' in kmer: continue
        counts[kmer] = counts.get(kmer, 0) + 1
        total += 1
    return counts, total

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    genomes = [
        "ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
        "pseudomonas_aeruginosa", "saureus_nctc8325",
        "mycobacterium_tuberculosis", "synechocystis_pcc6803",
        "bacteroides_fragilis", "borrelia_burgdorferi"
    ]

    # Collect all integrases and terminases
    all_integrases = []
    all_terminases = []
    all_other_phage = []

    for name in genomes:
        gff = f"{data_dir}/{name}.gff"
        fasta = f"{data_dir}/{name}.fasta"
        if not os.path.exists(gff) or not os.path.exists(fasta): continue

        integrases = extract_marker_genes(gff, fasta, ['integrase'])
        terminases = extract_marker_genes(gff, fasta, ['terminase'])
        other_phage = extract_marker_genes(gff, fasta, ['capsid', 'portal', 'tail', 'baseplate'])

        all_integrases.extend(integrases)
        all_terminases.extend(terminases)
        all_other_phage.extend(other_phage)

    print(f"Collected: {len(all_integrases)} integrases, {len(all_terminases)} terminases, {len(all_other_phage)} structural")

    # Length distribution
    int_lens = [g['length'] for g in all_integrases]
    term_lens = [g['length'] for g in all_terminases]

    import statistics
    if int_lens:
        print(f"\nIntegrase lengths: min={min(int_lens)}, median={statistics.median(int_lens):.0f}, max={max(int_lens)}")
    if term_lens:
        print(f"Terminase lengths: min={min(term_lens)}, median={statistics.median(term_lens):.0f}, max={max(term_lens)}")

    # Can we detect integrases by protein motifs?
    # Integrases have conserved catalytic residues: H-R-Y triad
    # But we're at DNA level. Let's look at amino acid patterns.
    print(f"\n=== Amino acid patterns in integrases ===")
    aa_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    def translate(dna):
        protein = []
        for i in range(0, len(dna) - 2, 3):
            codon = dna[i:i+3]
            protein.append(aa_table.get(codon, 'X'))
        return ''.join(protein)

    # Translate all integrases, find common short motifs (4-5 aa)
    proteins = [translate(g['seq']) for g in all_integrases if g['length'] >= 300]
    print(f"  {len(proteins)} integrase proteins (≥100aa)")

    # Find conserved 5-mer amino acid patterns
    motif_counts = {}
    for prot in proteins:
        seen = set()
        for i in range(len(prot) - 4):
            motif = prot[i:i+5]
            if '*' in motif or 'X' in motif: continue
            if motif not in seen:
                motif_counts[motif] = motif_counts.get(motif, 0) + 1
                seen.add(motif)

    # Most common motifs (present in >30% of integrases)
    threshold = len(proteins) * 0.3
    common = {k: v for k, v in motif_counts.items() if v >= threshold}
    common_sorted = sorted(common.items(), key=lambda x: -x[1])

    print(f"  5-mer motifs in ≥30% of integrases ({len(common_sorted)}):")
    for motif, count in common_sorted[:20]:
        pct = count * 100 // len(proteins)
        print(f"    {motif}: {count}/{len(proteins)} ({pct}%)")

    # Now check: are these motifs specific to integrases?
    # Check how many normal genes have the same motifs
    # Use dump-all-orfs TSV if available
    print(f"\n=== Motif specificity check ===")
    print(f"  (need to check against non-integrase genes)")

    # Save integrase DNA sequences for hex profile building
    all_int_dna = ''.join(g['seq'] for g in all_integrases if g['length'] >= 300)
    print(f"\n  Total integrase DNA: {len(all_int_dna)} bp from {len([g for g in all_integrases if g['length'] >= 300])} genes")
    print(f"  Enough for hex model: {'YES' if len(all_int_dna) > 10000 else 'NO'}")

    # Write integrase sequences to file for downstream use
    out_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'integrase_seqs.fasta')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        for i, g in enumerate(all_integrases):
            if g['length'] < 300: continue
            f.write(f">integrase_{i} {g['product'][:50]} len={g['length']}\n")
            seq = g['seq']
            for j in range(0, len(seq), 80):
                f.write(seq[j:j+80] + '\n')
    print(f"  Saved to {out_path}")

    # Also save terminases
    out_path2 = os.path.join(os.path.dirname(out_path), 'terminase_seqs.fasta')
    with open(out_path2, 'w') as f:
        for i, g in enumerate(all_terminases):
            if g['length'] < 300: continue
            f.write(f">terminase_{i} {g['product'][:50]} len={g['length']}\n")
            seq = g['seq']
            for j in range(0, len(seq), 80):
                f.write(seq[j:j+80] + '\n')
    print(f"  Saved terminases to {out_path2}")

if __name__ == "__main__":
    main()
