#!/usr/bin/env python3
"""Leave-one-out validation of anchor caller.
For each genome: build library WITHOUT it, then find anchors IN it."""

import os, json, statistics
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

def extract_all_proteins(gff_path, fasta_path):
    """Extract ALL CDS proteins with gene names."""
    seqs = read_fasta(fasta_path)
    proteins = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS': continue
            if 'pseudo=true' in fields[8]: continue
            seqid = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            name = None
            for part in fields[8].split(';'):
                if part.startswith('gene='): name = part.split('=')[1]
            if seqid not in seqs: continue
            dna = seqs[seqid][start:end]
            if strand == '-': dna = revcomp(dna)
            prot = translate(dna)
            if len(prot) < 30: continue
            proteins.append({
                'gene_name': name,
                'start': start + 1, 'end': end, 'strand': strand,
                'protein': prot, 'length': len(prot),
            })
    return proteins

def build_family_fingerprint(proteins, k=5):
    """Build fingerprint: k-mers present in ≥50% of family members."""
    kmer_counts = defaultdict(int)
    for prot in proteins:
        seen = set()
        for i in range(len(prot) - k + 1):
            kmer = prot[i:i+k]
            if 'X' not in kmer and kmer not in seen:
                kmer_counts[kmer] += 1
                seen.add(kmer)
    threshold = max(len(proteins) * 0.4, 2)  # present in ≥40% or ≥2 members
    fingerprint = set(km for km, cnt in kmer_counts.items() if cnt >= threshold)
    return fingerprint

def score_protein_fingerprint(prot, fingerprint, k=5):
    """Score = fraction of protein k-mers that match fingerprint."""
    hits = 0
    total = 0
    for i in range(len(prot) - k + 1):
        kmer = prot[i:i+k]
        if 'X' not in kmer:
            total += 1
            if kmer in fingerprint:
                hits += 1
    return hits / max(total, 1)

def build_library_without(anchor_families, exclude_genome, all_proteins_by_genome, data_dir):
    """Build fingerprints excluding one genome."""
    profiles = {}
    for family_name, family_genomes in anchor_families.items():
        family_prots = []
        for gname in family_genomes:
            if gname == exclude_genome: continue
            for p in all_proteins_by_genome.get(gname, []):
                if p['gene_name'] == family_name:
                    family_prots.append(p['protein'])
                    break
        if len(family_prots) < 3: continue
        fingerprint = build_family_fingerprint(family_prots)
        if len(fingerprint) < 5: continue  # too few specific k-mers
        median_len = statistics.median([len(p) for p in family_prots])
        profiles[family_name] = {
            'profile': fingerprint,
            'median_length': median_len,
            'n_training': len(family_prots),
        }
    return profiles

def call_anchors(proteins, profiles, length_tolerance=0.5):
    """Find best anchor assignment for each protein using fingerprint matching."""
    assignments = {}

    for family_name, prof_info in profiles.items():
        fingerprint = prof_info['profile']
        med_len = prof_info['median_length']

        best_score = -1
        best_prot = None
        second_score = -1

        for p in proteins:
            if abs(p['length'] - med_len) > med_len * length_tolerance:
                continue
            score = score_protein_fingerprint(p['protein'], fingerprint)
            if score > best_score:
                second_score = best_score
                best_score = score
                best_prot = p
            elif score > second_score:
                second_score = score

        if best_prot is not None and best_score > 0.02:  # at least 2% k-mers match
            specificity = best_score / max(second_score, 1e-10) if second_score > 0 else 100.0
            assignments[family_name] = {
                'protein': best_prot,
                'score': best_score,
                'specificity': specificity,
            }

    return assignments

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"
    lib_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'anchor_library')

    genomes = [
        "ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
        "pseudomonas_aeruginosa", "saureus_nctc8325",
        "mycobacterium_tuberculosis", "mycoplasma_genitalium",
        "synechocystis_pcc6803", "bacteroides_fragilis",
        "borrelia_burgdorferi"
    ]

    # Load anchor library metadata
    with open(os.path.join(lib_dir, 'anchor_library.json')) as f:
        lib_meta = json.load(f)

    # Load anchor graph to get which genomes each anchor is in
    graph_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'anchor_graph.json')
    with open(graph_path) as f:
        graph = json.load(f)

    anchor_families = {}  # family_name -> list of genomes
    for name in lib_meta:
        if name in graph['anchors']:
            anchor_families[name] = graph['anchors'][name]

    # Pre-load all proteins
    print("Loading proteins...")
    all_proteins = {}
    for gname in genomes:
        gff = f"{data_dir}/{gname}.gff"
        fasta = f"{data_dir}/{gname}.fasta"
        if os.path.exists(gff):
            prots = extract_all_proteins(gff, fasta)
            all_proteins[gname] = prots
            print(f"  {gname}: {len(prots)} proteins")

    # Leave-one-out validation
    print(f"\n{'='*70}")
    print(f"  LEAVE-ONE-OUT VALIDATION")
    print(f"{'='*70}")
    print(f"\n{'Genome':<25} {'Expected':>8} {'Found':>5} {'Correct':>7} {'Wrong':>5} {'Missed':>6} {'Recall':>7} {'Precision':>9}")
    print("-" * 80)

    total_expected = 0
    total_correct = 0
    total_wrong = 0
    total_missed = 0

    for exclude in genomes:
        # Build library without this genome
        profiles = build_library_without(anchor_families, exclude, all_proteins, data_dir)

        # How many anchors should this genome have?
        expected_anchors = set()
        for name, glist in anchor_families.items():
            if exclude in glist and name in profiles:
                expected_anchors.add(name)

        # Call anchors in excluded genome
        target_prots = all_proteins.get(exclude, [])
        assignments = call_anchors(target_prots, profiles)

        # Evaluate
        correct = 0
        wrong = 0
        for family_name, assignment in assignments.items():
            called_gene = assignment['protein']['gene_name']
            if called_gene == family_name:
                correct += 1
            else:
                wrong += 1

        found = len(assignments)
        missed = len(expected_anchors) - correct
        recall = correct / max(len(expected_anchors), 1)
        precision = correct / max(found, 1)

        print(f"{exclude:<25} {len(expected_anchors):>8} {found:>5} {correct:>7} {wrong:>5} {missed:>6} {recall:>7.1%} {precision:>9.1%}")

        total_expected += len(expected_anchors)
        total_correct += correct
        total_wrong += wrong
        total_missed += missed

    print("-" * 80)
    overall_recall = total_correct / max(total_expected, 1)
    overall_precision = total_correct / max(total_correct + total_wrong, 1)
    print(f"{'TOTAL':<25} {total_expected:>8} {total_correct+total_wrong:>5} {total_correct:>7} {total_wrong:>5} {total_missed:>6} {overall_recall:>7.1%} {overall_precision:>9.1%}")

if __name__ == "__main__":
    main()
