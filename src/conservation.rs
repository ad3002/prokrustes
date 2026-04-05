//! Cross-genome conservation scoring (de novo comparative genomics).
//!
//! Compares predicted proteins across multiple genomes using protein k-mer
//! similarity. No external databases — only our own predictions.
//!
//! Each gene gets a conservation level:
//!   "conserved" — ortholog in ≥1 comparison genome (FP rate ~1%)
//!   "unique"    — no ortholog found (FP rate ~15%)

use std::collections::HashMap;
use crate::io;
use crate::types::Gene;

const K: usize = 9; // protein k-mer size (9 with Dayhoff-6: 6^9=10M variants)

/// Dayhoff-6 reduced alphabet: groups conservative substitutions.
/// C=0, AGPST=1, DENQ=2, HKR=3, ILMV=4, FWY=5
/// This captures I↔V, D↔E, K↔R etc. as same "letter".
#[inline]
fn dayhoff6(aa: u8) -> u8 {
    match aa {
        b'C' => 0,
        b'A' | b'G' | b'P' | b'S' | b'T' => 1,
        b'D' | b'E' | b'N' | b'Q' => 2,
        b'H' | b'K' | b'R' => 3,
        b'I' | b'L' | b'M' | b'V' => 4,
        b'F' | b'W' | b'Y' => 5,
        _ => 7, // unknown
    }
}

/// Convert protein to reduced alphabet
fn reduce(protein: &[u8]) -> Vec<u8> {
    protein.iter().map(|&aa| dayhoff6(aa)).collect()
}

static CODON_TABLE: &[(u8, u8, u8, u8)] = &[
    (b'T',b'T',b'T',b'F'),(b'T',b'T',b'C',b'F'),(b'T',b'T',b'A',b'L'),(b'T',b'T',b'G',b'L'),
    (b'C',b'T',b'T',b'L'),(b'C',b'T',b'C',b'L'),(b'C',b'T',b'A',b'L'),(b'C',b'T',b'G',b'L'),
    (b'A',b'T',b'T',b'I'),(b'A',b'T',b'C',b'I'),(b'A',b'T',b'A',b'I'),(b'A',b'T',b'G',b'M'),
    (b'G',b'T',b'T',b'V'),(b'G',b'T',b'C',b'V'),(b'G',b'T',b'A',b'V'),(b'G',b'T',b'G',b'V'),
    (b'T',b'C',b'T',b'S'),(b'T',b'C',b'C',b'S'),(b'T',b'C',b'A',b'S'),(b'T',b'C',b'G',b'S'),
    (b'C',b'C',b'T',b'P'),(b'C',b'C',b'C',b'P'),(b'C',b'C',b'A',b'P'),(b'C',b'C',b'G',b'P'),
    (b'A',b'C',b'T',b'T'),(b'A',b'C',b'C',b'T'),(b'A',b'C',b'A',b'T'),(b'A',b'C',b'G',b'T'),
    (b'G',b'C',b'T',b'A'),(b'G',b'C',b'C',b'A'),(b'G',b'C',b'A',b'A'),(b'G',b'C',b'G',b'A'),
    (b'T',b'A',b'T',b'Y'),(b'T',b'A',b'C',b'Y'),
    (b'C',b'A',b'T',b'H'),(b'C',b'A',b'C',b'H'),(b'C',b'A',b'A',b'Q'),(b'C',b'A',b'G',b'Q'),
    (b'A',b'A',b'T',b'N'),(b'A',b'A',b'C',b'N'),(b'A',b'A',b'A',b'K'),(b'A',b'A',b'G',b'K'),
    (b'G',b'A',b'T',b'D'),(b'G',b'A',b'C',b'D'),(b'G',b'A',b'A',b'E'),(b'G',b'A',b'G',b'E'),
    (b'T',b'G',b'T',b'C'),(b'T',b'G',b'C',b'C'),(b'T',b'G',b'G',b'W'),
    (b'C',b'G',b'T',b'R'),(b'C',b'G',b'C',b'R'),(b'C',b'G',b'A',b'R'),(b'C',b'G',b'G',b'R'),
    (b'A',b'G',b'T',b'S'),(b'A',b'G',b'C',b'S'),(b'A',b'G',b'A',b'R'),(b'A',b'G',b'G',b'R'),
    (b'G',b'G',b'T',b'G'),(b'G',b'G',b'C',b'G'),(b'G',b'G',b'A',b'G'),(b'G',b'G',b'G',b'G'),
];

fn translate_codon(a: u8, b: u8, c: u8) -> Option<u8> {
    for &(x, y, z, aa) in CODON_TABLE {
        if a == x && b == y && c == z { return Some(aa); }
    }
    None
}

fn translate(dna: &[u8]) -> Vec<u8> {
    let mut prot = Vec::new();
    let mut i = 0;
    while i + 3 <= dna.len() {
        if let Some(aa) = translate_codon(dna[i], dna[i+1], dna[i+2]) {
            prot.push(aa);
        } else {
            break;
        }
        i += 3;
    }
    prot
}

fn get_protein(genome: &[u8], gene: &Gene) -> Vec<u8> {
    let seq = if gene.is_plus {
        genome[gene.start.saturating_sub(1)..gene.end.min(genome.len())].to_vec()
    } else {
        io::rev_comp(&genome[gene.start.saturating_sub(1)..gene.end.min(genome.len())])
    };
    translate(&seq)
}

/// K-mer index for fast protein comparison.
/// Uses Dayhoff-6 reduced alphabet for substitution-tolerant matching.
pub struct KmerIndex {
    index: HashMap<[u8; K], Vec<u16>>, // reduced kmer → protein IDs
    n_proteins: usize,
}

impl KmerIndex {
    /// Build index from predicted genes on a comparison genome.
    pub fn build(genome: &[u8], genes: &[Gene]) -> Self {
        let mut index: HashMap<[u8; K], Vec<u16>> = HashMap::new();
        let n = genes.len().min(u16::MAX as usize);
        for (i, gene) in genes.iter().enumerate().take(n) {
            let prot = get_protein(genome, gene);
            let reduced = reduce(&prot);
            if reduced.len() < K { continue; }
            for j in 0..reduced.len() - K + 1 {
                let mut kmer = [0u8; K];
                kmer.copy_from_slice(&reduced[j..j+K]);
                if kmer.contains(&7) { continue; } // skip unknown
                index.entry(kmer).or_default().push(i as u16);
            }
        }
        KmerIndex { index, n_proteins: n }
    }

    /// Score a protein against this index (using reduced alphabet).
    /// Returns fraction of k-mers shared with best hit (0.0 to 1.0).
    pub fn score(&self, protein: &[u8]) -> f64 {
        let reduced = reduce(protein);
        if reduced.len() < K { return 0.0; }
        let total = reduced.len() - K + 1;
        let mut hits = vec![0u16; self.n_proteins];
        for j in 0..total {
            let mut kmer = [0u8; K];
            kmer.copy_from_slice(&reduced[j..j+K]);
            if kmer.contains(&7) { continue; }
            if let Some(ids) = self.index.get(&kmer) {
                for &id in ids {
                    hits[id as usize] += 1;
                }
            }
        }
        let best = hits.iter().max().copied().unwrap_or(0) as f64;
        best / total as f64
    }
}

/// Predict genes on a comparison genome and build k-mer index.
pub fn build_comparison_index(fasta_path: &str) -> Option<(Vec<Gene>, KmerIndex)> {
    let genome = io::read_fasta(fasta_path);
    if genome.len() < 1000 { return None; }
    let (genes, _) = crate::pipeline::annotate(&genome);
    let index = KmerIndex::build(&genome, &genes);
    eprintln!("  Comparison {}: {} genes, {} k-mers",
        fasta_path.rsplit('/').next().unwrap_or(fasta_path),
        genes.len(), index.index.len());
    Some((genes, index))
}

/// Score conservation of target genes against multiple comparison indexes.
/// Returns conservation score for each gene (0.0 = unique, 1.0 = highly conserved).
pub fn score_conservation(
    genome: &[u8],
    genes: &[Gene],
    indexes: &[&KmerIndex],
) -> Vec<f64> {
    genes.iter().map(|gene| {
        let prot = get_protein(genome, gene);
        if prot.len() < K { return 0.0; }
        // Best conservation across all comparison genomes
        let mut best = 0.0f64;
        for idx in indexes {
            let s = idx.score(&prot);
            if s > best { best = s; }
        }
        best
    }).collect()
}

/// Classify conservation level.
pub fn conservation_label(score: f64) -> &'static str {
    if score > 0.3 { "conserved" }
    else if score > 0.1 { "weak" }
    else { "unique" }
}

/// Find best ortholog length for a gene.
/// Returns (conservation_score, ortholog_length) or None.
pub fn find_ortholog_length(protein: &[u8], indexes: &[&KmerIndex], comp_genomes: &[(&[u8], &[Gene])]) -> Option<(f64, usize)> {
    let reduced = reduce(protein);
    if reduced.len() < K { return None; }
    let total = reduced.len() - K + 1;

    for (idx_i, idx) in indexes.iter().enumerate() {
        let mut hits = vec![0u16; idx.n_proteins];
        for j in 0..total {
            let mut kmer = [0u8; K];
            kmer.copy_from_slice(&reduced[j..j+K]);
            if kmer.contains(&7) { continue; }
            if let Some(ids) = idx.index.get(&kmer) {
                for &id in ids { hits[id as usize] += 1; }
            }
        }
        let best_id = hits.iter().enumerate().max_by_key(|(_,&c)| c).map(|(i,_)| i);
        if let Some(bi) = best_id {
            let best_frac = hits[bi] as f64 / total as f64;
            if best_frac > 0.3 {
                if idx_i < comp_genomes.len() {
                    let (_, comp_genes) = comp_genomes[idx_i];
                    if bi < comp_genes.len() {
                        return Some((best_frac, comp_genes[bi].length));
                    }
                }
            }
        }
    }
    None
}

/// Conservation-guided start correction.
/// If our gene is significantly shorter than its ortholog, try extending
/// upstream to find a start codon that gives a length closer to the ortholog.
pub fn correct_starts_by_orthologs(
    genome: &[u8],
    genes: &mut [Gene],
    all_orfs: &[Gene],
    indexes: &[&KmerIndex],
    comp_genomes: &[(&[u8], &[Gene])],
) -> usize {
    let mut corrected = 0;

    for gene in genes.iter_mut() {
        let prot = get_protein(genome, gene);
        let orth = find_ortholog_length(&prot, indexes, comp_genomes);
        if orth.is_none() { continue; }
        let (_, orth_len) = orth.unwrap();

        // Only correct if ortholog is LONGER (>10% or >90bp)
        if orth_len <= gene.length { continue; }
        let diff = orth_len - gene.length;
        if diff < 90 && diff < gene.length / 10 { continue; }

        let target_len = orth_len;
        if diff > 2000 { continue; }

        // Find a longer ORF in the same stop group that's closer to ortholog length
        let mut best_alt: Option<&Gene> = None;
        let mut best_diff = diff;

        for orf in all_orfs.iter() {
            if orf.is_plus != gene.is_plus { continue; }
            if orf.end != gene.end { continue; } // same stop codon
            if orf.length <= gene.length { continue; } // must be longer
            if orf.length < 90 { continue; }

            let new_diff = if orf.length > target_len {
                orf.length - target_len
            } else {
                target_len - orf.length
            };

            if new_diff < best_diff {
                best_diff = new_diff;
                best_alt = Some(orf);
            }
        }

        // Accept correction if new length is within 20% of ortholog
        if let Some(alt) = best_alt {
            if best_diff < target_len / 5 {
                gene.start = alt.start;
                gene.end = alt.end;
                gene.length = alt.length;
                gene.seq_start = alt.seq_start;
                gene.seq_end = alt.seq_end;
                gene.start_codon = alt.start_codon;
                corrected += 1;
            }
        }
    }

    corrected
}
