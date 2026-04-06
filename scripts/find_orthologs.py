#!/usr/bin/env python3
"""Find single-copy orthologous genes across 10 genomes.
Method: translate all CDS → protein sequences → find shared protein k-mers
→ cluster by reciprocal best k-mer hits → count single-copy orthologs."""

import os, sys
from collections import defaultdict

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
    return ''.join(aa_table.get(dna[i:i+3], 'X') for i in range(0, len(dna) - 2, 3)).split('*')[0]

def read_fasta(path):
    seqs = {}
    name = ""
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name: seqs[name] = ''.join(parts)
                name = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip().upper())
    if name: seqs[name] = ''.join(parts)
    return seqs

def revcomp(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp.get(c, 'N') for c in reversed(seq))

def extract_proteins(gff_path, fasta_path):
    """Extract protein sequences from CDS annotations."""
    seqs = read_fasta(fasta_path)
    proteins = []

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS': continue
            attrs = fields[8]
            if 'pseudo=true' in attrs or 'pseudogene=' in attrs: continue

            seqid = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]

            if seqid not in seqs: continue
            dna = seqs[seqid][start:end]
            if strand == '-':
                dna = revcomp(dna)

            prot = translate(dna)
            if len(prot) < 50: continue  # skip short

            # Get gene name/locus_tag
            gene_id = ""
            for part in attrs.split(';'):
                if part.startswith('locus_tag='):
                    gene_id = part.split('=')[1]
                    break
                if part.startswith('Name=') and not gene_id:
                    gene_id = part.split('=')[1]

            proteins.append({
                'id': gene_id or f"{seqid}:{start}-{end}",
                'start': start, 'end': end, 'strand': strand,
                'protein': prot, 'length': len(prot),
            })

    return proteins

def protein_kmers(prot, k=5):
    """Extract amino acid k-mers from protein sequence."""
    kmers = set()
    for i in range(len(prot) - k + 1):
        kmer = prot[i:i+k]
        if 'X' not in kmer and '*' not in kmer:
            kmers.add(kmer)
    return kmers

def jaccard(set1, set2):
    if not set1 or not set2: return 0.0
    return len(set1 & set2) / len(set1 | set2)

def find_best_hits(proteins_a, proteins_b, k=5, min_jaccard=0.15):
    """Find reciprocal best hits by protein k-mer Jaccard similarity."""
    # Pre-compute k-mers
    kmers_a = [(p, protein_kmers(p['protein'], k)) for p in proteins_a]
    kmers_b = [(p, protein_kmers(p['protein'], k)) for p in proteins_b]

    # Best hit A → B
    best_ab = {}
    for pa, ka in kmers_a:
        best_j = 0
        best_p = None
        for pb, kb in kmers_b:
            # Quick length filter
            if abs(pa['length'] - pb['length']) > max(pa['length'], pb['length']) * 0.5:
                continue
            j = jaccard(ka, kb)
            if j > best_j:
                best_j = j
                best_p = pb
        if best_j >= min_jaccard:
            best_ab[pa['id']] = (best_p['id'], best_j)

    # Best hit B → A
    best_ba = {}
    for pb, kb in kmers_b:
        best_j = 0
        best_p = None
        for pa, ka in kmers_a:
            if abs(pa['length'] - pb['length']) > max(pa['length'], pb['length']) * 0.5:
                continue
            j = jaccard(ka, kb)
            if j > best_j:
                best_j = j
                best_p = pa
        if best_j >= min_jaccard:
            best_ba[pb['id']] = (best_p['id'], best_j)

    # Reciprocal best hits
    rbh = []
    for a_id, (b_id, j_ab) in best_ab.items():
        if b_id in best_ba and best_ba[b_id][0] == a_id:
            rbh.append((a_id, b_id, j_ab))

    return rbh

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"

    genomes = [
        ("ecoli_k12", "Escherichia coli"),
        ("salmonella_lt2", "Salmonella LT2"),
        ("bacillus_subtilis", "Bacillus subtilis"),
        ("pseudomonas_aeruginosa", "Pseudomonas aeruginosa"),
        ("saureus_nctc8325", "S. aureus"),
        ("mycobacterium_tuberculosis", "M. tuberculosis"),
        ("mycoplasma_genitalium", "Mycoplasma genitalium"),
        ("synechocystis_pcc6803", "Synechocystis PCC 6803"),
        ("bacteroides_fragilis", "Bacteroides fragilis"),
        ("borrelia_burgdorferi", "Borrelia burgdorferi"),
    ]

    # Extract proteins from each genome
    all_proteins = {}
    for name, label in genomes:
        gff = f"{data_dir}/{name}.gff"
        fasta = f"{data_dir}/{name}.fasta"
        if not os.path.exists(gff): continue
        prots = extract_proteins(gff, fasta)
        all_proteins[name] = prots
        print(f"{label}: {len(prots)} proteins (≥50aa)")

    print(f"\nTotal: {sum(len(v) for v in all_proteins.values())} proteins")

    # Find pairwise orthologs between close species first
    # (full all-vs-all is O(n²) and slow)
    pairs_to_check = [
        ("ecoli_k12", "salmonella_lt2"),         # Enterobacteria
        ("ecoli_k12", "bacillus_subtilis"),       # Firmicutes vs Gamma
        ("ecoli_k12", "pseudomonas_aeruginosa"),  # Gamma-proteobacteria
        ("ecoli_k12", "saureus_nctc8325"),        # Firmicute
        ("ecoli_k12", "mycobacterium_tuberculosis"), # Actinobacteria
        ("ecoli_k12", "mycoplasma_genitalium"),   # Tenericutes
        ("ecoli_k12", "synechocystis_pcc6803"),   # Cyanobacteria
        ("ecoli_k12", "bacteroides_fragilis"),     # Bacteroidetes
        ("ecoli_k12", "borrelia_burgdorferi"),     # Spirochaetes
    ]

    print(f"\n{'='*60}")
    print(f"  PAIRWISE ORTHOLOG SEARCH")
    print(f"{'='*60}")

    pair_orthologs = {}
    for a_name, b_name in pairs_to_check:
        if a_name not in all_proteins or b_name not in all_proteins:
            continue

        print(f"\n{a_name} vs {b_name}...", end=" ", flush=True)
        rbh = find_best_hits(all_proteins[a_name], all_proteins[b_name])
        pair_orthologs[(a_name, b_name)] = rbh

        if rbh:
            jaccards = [j for _, _, j in rbh]
            import statistics
            print(f"{len(rbh)} orthologs (median similarity: {statistics.median(jaccards):.3f})")
        else:
            print("0 orthologs")

    # Find genes shared across multiple genomes
    print(f"\n{'='*60}")
    print(f"  MULTI-GENOME ORTHOLOGS")
    print(f"{'='*60}")

    # Count how many genomes each E. coli gene has orthologs in
    ecoli_ortho_count = defaultdict(int)
    ecoli_ortho_genomes = defaultdict(set)

    for (a, b), rbhs in pair_orthologs.items():
        for ecoli_id, other_id, sim in rbhs:
            ecoli_ortho_count[ecoli_id] += 1
            ecoli_ortho_genomes[ecoli_id].add(b)

    # Distribution of conservation
    counts = list(ecoli_ortho_count.values())
    for n_genomes in range(1, 10):
        n_genes = sum(1 for c in counts if c >= n_genomes)
        print(f"  E. coli genes with orthologs in ≥{n_genomes} genomes: {n_genes}")

    # Core genes (present in all 9 comparisons)
    core = [gid for gid, c in ecoli_ortho_count.items() if c >= 8]
    print(f"\n  CORE GENES (orthologs in ≥8/9 genomes): {len(core)}")

    # Show some examples
    if core:
        for gid in core[:10]:
            n = ecoli_ortho_count[gid]
            genomes_str = ', '.join(sorted(ecoli_ortho_genomes[gid]))
            print(f"    {gid}: {n} orthologs")

if __name__ == "__main__":
    main()
