#!/usr/bin/env python3
"""Build reusable anchor library: select best families, extract proteins,
compute MSA quality metrics, build simple profile models."""

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

def extract_protein_by_gene_name(gff_path, fasta_path, target_name):
    """Extract protein sequence for a specific gene name."""
    seqs = read_fasta(fasta_path)
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS': continue
            if 'pseudo=true' in fields[8]: continue
            if f'gene={target_name}' not in fields[8] and f';gene={target_name};' not in fields[8]:
                # More careful check
                has_gene = False
                for part in fields[8].split(';'):
                    if part == f'gene={target_name}':
                        has_gene = True
                        break
                if not has_gene: continue

            seqid = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            if seqid not in seqs: continue
            dna = seqs[seqid][start:end]
            if strand == '-': dna = revcomp(dna)
            prot = translate(dna)
            if len(prot) >= 30:
                return prot
    return None

def protein_identity(p1, p2):
    """Simple identity: shared k-mers / union k-mers."""
    k = 3
    s1 = set(p1[i:i+k] for i in range(len(p1)-k+1))
    s2 = set(p2[i:i+k] for i in range(len(p2)-k+1))
    if not s1 or not s2: return 0
    return len(s1 & s2) / len(s1 | s2)

def build_kmer_profile(proteins, k=3):
    """Build amino acid k-mer frequency profile from multiple proteins."""
    counts = defaultdict(float)
    total = 0
    for prot in proteins:
        for i in range(len(prot) - k + 1):
            kmer = prot[i:i+k]
            if 'X' not in kmer:
                counts[kmer] += 1
                total += 1
    # Normalize
    profile = {}
    for kmer, cnt in counts.items():
        profile[kmer] = cnt / max(total, 1)
    return profile

def score_protein_against_profile(prot, profile, k=3):
    """Score a protein against a k-mer profile."""
    total = 0.0
    n = 0
    for i in range(len(prot) - k + 1):
        kmer = prot[i:i+k]
        if kmer in profile:
            total += profile[kmer]
            n += 1
    return total / max(n, 1)

def main():
    data_dir = "/Users/akomissarov/Dropbox/workspace/new/renorm/data"
    norm_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'normalized')
    out_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'anchor_library')
    os.makedirs(out_dir, exist_ok=True)

    # Load anchor graph
    graph_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'anchor_graph.json')
    with open(graph_path) as f:
        graph = json.load(f)

    genomes = [
        "ecoli_k12", "salmonella_lt2", "bacillus_subtilis",
        "pseudomonas_aeruginosa", "saureus_nctc8325",
        "mycobacterium_tuberculosis", "mycoplasma_genitalium",
        "synechocystis_pcc6803", "bacteroides_fragilis",
        "borrelia_burgdorferi"
    ]

    # Step 1: Select best anchor families
    # Criteria: ≥6 genomes, single-copy (check), consistent length
    print("=== SELECTING BEST ANCHOR FAMILIES ===\n")

    anchor_genomes = graph['anchors']  # name -> [genome_list]
    candidates = []

    for name, genome_list in anchor_genomes.items():
        n_genomes = len(genome_list)
        if n_genomes < 6: continue

        # Check single-copy: does this gene appear once per genome?
        single_copy = True
        for gname in genome_list:
            chain = graph['chains'].get(gname, [])
            copies = sum(1 for g in chain if g['name'] == name)
            if copies > 1:
                single_copy = False
                break

        if not single_copy: continue

        # Get protein lengths across genomes
        lengths = []
        for gname in genome_list:
            chain = graph['chains'].get(gname, [])
            for g in chain:
                if g['name'] == name:
                    lengths.append(g['length'] // 3)  # approx AA length

        if not lengths: continue
        med_len = statistics.median(lengths)
        len_cv = statistics.stdev(lengths) / max(med_len, 1) if len(lengths) > 1 else 0

        candidates.append({
            'name': name,
            'n_genomes': n_genomes,
            'median_length_aa': med_len,
            'length_cv': len_cv,
            'single_copy': True,
        })

    # Sort by: n_genomes desc, then length_cv asc
    candidates.sort(key=lambda c: (-c['n_genomes'], c['length_cv']))

    print(f"Candidates (≥6 genomes, single-copy): {len(candidates)}")
    print(f"\nTop 20:")
    for c in candidates[:20]:
        print(f"  {c['name']:>12}: {c['n_genomes']} genomes, {c['median_length_aa']:.0f}aa, CV={c['length_cv']:.2f}")

    # Select top 100 (or all if fewer)
    selected = candidates[:min(100, len(candidates))]
    print(f"\nSelected {len(selected)} anchors for library")

    # Step 2: Extract proteins and build profiles
    print(f"\n=== EXTRACTING PROTEINS ===\n")

    library = {}
    n_extracted = 0

    for anchor in selected:
        name = anchor['name']
        proteins = []

        for gname in genomes:
            # Use original (not normalized) FASTA + GFF for protein extraction
            gff = f"{data_dir}/{gname}.gff"
            fasta = f"{data_dir}/{gname}.fasta"
            if not os.path.exists(gff): continue

            prot = extract_protein_by_gene_name(gff, fasta, name)
            if prot:
                proteins.append({'genome': gname, 'sequence': prot, 'length': len(prot)})

        if len(proteins) < 4: continue

        # Build k-mer profile
        seqs = [p['sequence'] for p in proteins]
        profile = build_kmer_profile(seqs, k=4)

        # Quality: average pairwise identity
        identities = []
        for i in range(len(seqs)):
            for j in range(i+1, min(i+5, len(seqs))):
                identities.append(protein_identity(seqs[i], seqs[j]))

        avg_identity = statistics.mean(identities) if identities else 0

        library[name] = {
            'n_genomes': len(proteins),
            'median_length': statistics.median([p['length'] for p in proteins]),
            'avg_identity': avg_identity,
            'profile': profile,
            'proteins': proteins,
        }
        n_extracted += 1

    print(f"Built profiles for {n_extracted} anchor families")

    # Step 3: Quality statistics
    print(f"\n=== LIBRARY QUALITY ===\n")

    identities = [v['avg_identity'] for v in library.values()]
    lengths = [v['median_length'] for v in library.values()]
    n_genomes_list = [v['n_genomes'] for v in library.values()]

    print(f"Families: {len(library)}")
    print(f"Identity: median={statistics.median(identities):.3f}, min={min(identities):.3f}")
    print(f"Length:   median={statistics.median(lengths):.0f}aa")
    print(f"Coverage: median={statistics.median(n_genomes_list)} genomes")

    # By identity buckets
    for thresh in [0.1, 0.2, 0.3, 0.5]:
        n = sum(1 for x in identities if x >= thresh)
        print(f"  Identity ≥{thresh}: {n}/{len(identities)}")

    # Save library
    # Save profiles as JSON (without full protein sequences for size)
    lib_json = {}
    for name, info in library.items():
        lib_json[name] = {
            'n_genomes': info['n_genomes'],
            'median_length': info['median_length'],
            'avg_identity': info['avg_identity'],
            'profile_size': len(info['profile']),
        }

    lib_path = os.path.join(out_dir, 'anchor_library.json')
    with open(lib_path, 'w') as f:
        json.dump(lib_json, f, indent=2)
    print(f"\nLibrary metadata saved: {lib_path}")

    # Save protein sequences as FASTA per family
    for name, info in library.items():
        fasta_path = os.path.join(out_dir, f"{name}.fasta")
        with open(fasta_path, 'w') as f:
            for p in info['proteins']:
                f.write(f">{name}_{p['genome']} len={p['length']}\n")
                seq = p['sequence']
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')

    print(f"Protein FASTA files: {out_dir}/ ({len(library)} families)")

    # Save k-mer profiles
    profiles_path = os.path.join(out_dir, 'profiles.json')
    profiles_json = {}
    for name, info in library.items():
        profiles_json[name] = {
            'profile': info['profile'],
            'n_genomes': info['n_genomes'],
            'median_length': info['median_length'],
        }
    with open(profiles_path, 'w') as f:
        json.dump(profiles_json, f)
    print(f"K-mer profiles saved: {profiles_path}")

if __name__ == "__main__":
    main()
