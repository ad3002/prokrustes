//! Learned start codon context model.
//!
//! Takes a fixed window around a candidate start codon, encodes nucleotides
//! as one-hot features, and classifies: "is this a real translation start?"
//!
//! Self-trained on confident predictions from the same genome.
//! Taxon-specific (Enterobacteria share patterns) but not species-specific.

use crate::types::Gene;
use crate::io::nt4;

/// Window: 60bp upstream + start codon (3bp) + 60bp downstream = 123bp total
/// Upstream captures: SD motif (5-15bp), promoter (-35bp), AT context
/// Downstream captures: first 20 codons (codon ramp, initial coding pattern)
const UP: usize = 60;
const DOWN: usize = 60;
const WIN: usize = UP + 3 + DOWN; // 123
const N_FEAT: usize = WIN * 4;    // 492 features (one-hot)

/// Logistic regression model for start context.
pub struct StartModel {
    weights: Vec<f64>, // N_FEAT weights
    bias: f64,
}

impl StartModel {
    /// Extract one-hot feature vector for a candidate start position.
    /// Returns None if window extends past sequence boundaries or has ambiguous bases.
    fn extract(seq: &[u8], start_pos: usize) -> Option<Vec<f64>> {
        if start_pos < UP { return None; }
        let end = start_pos + 3 + DOWN;
        if end > seq.len() { return None; }

        let window = &seq[start_pos - UP..end];
        let mut feat = vec![0.0f64; N_FEAT];

        for (i, &b) in window.iter().enumerate() {
            let n = nt4(b);
            if n > 3 { return None; } // ambiguous base
            feat[i * 4 + n] = 1.0;
        }
        Some(feat)
    }

    /// Sigmoid function.
    fn sigmoid(x: f64) -> f64 {
        1.0 / (1.0 + (-x).exp())
    }

    /// Score a candidate start position.
    pub fn score(&self, seq: &[u8], start_pos: usize) -> f64 {
        match Self::extract(seq, start_pos) {
            Some(feat) => {
                let z: f64 = feat.iter().zip(self.weights.iter())
                    .map(|(f, w)| f * w)
                    .sum::<f64>() + self.bias;
                Self::sigmoid(z)
            }
            None => 0.5, // neutral if can't extract
        }
    }

    /// Train model on positive (real starts) and negative (wrong starts) examples.
    /// Uses SGD with L2 regularization.
    pub fn train(
        seq: &[u8],
        positives: &[usize], // seq_start positions of real starts
        negatives: &[usize], // seq_start positions of non-starts
        n_epochs: usize,
        lr: f64,
        l2: f64,
    ) -> Option<StartModel> {
        if positives.len() < 50 || negatives.len() < 50 { return None; }

        // Collect valid training examples
        let mut pos_feats: Vec<Vec<f64>> = Vec::new();
        let mut neg_feats: Vec<Vec<f64>> = Vec::new();

        for &p in positives {
            if let Some(f) = Self::extract(seq, p) {
                pos_feats.push(f);
            }
        }
        for &p in negatives {
            if let Some(f) = Self::extract(seq, p) {
                neg_feats.push(f);
            }
        }

        if pos_feats.len() < 40 || neg_feats.len() < 40 { return None; }

        let mut weights = vec![0.0f64; N_FEAT];
        let mut bias = 0.0f64;

        // Balance classes: subsample the larger class
        let n_per_class = pos_feats.len().min(neg_feats.len());

        for _epoch in 0..n_epochs {
            // Positive examples (label = 1)
            for feat in pos_feats.iter().take(n_per_class) {
                let z: f64 = feat.iter().zip(weights.iter())
                    .map(|(f, w)| f * w).sum::<f64>() + bias;
                let pred = Self::sigmoid(z);
                let err = pred - 1.0; // gradient for label=1

                for (j, &f) in feat.iter().enumerate() {
                    weights[j] -= lr * (err * f + l2 * weights[j]);
                }
                bias -= lr * err;
            }

            // Negative examples (label = 0)
            for feat in neg_feats.iter().take(n_per_class) {
                let z: f64 = feat.iter().zip(weights.iter())
                    .map(|(f, w)| f * w).sum::<f64>() + bias;
                let pred = Self::sigmoid(z);
                let err = pred - 0.0; // gradient for label=0

                for (j, &f) in feat.iter().enumerate() {
                    weights[j] -= lr * (err * f + l2 * weights[j]);
                }
                bias -= lr * err;
            }
        }

        Some(StartModel { weights, bias })
    }
}

/// Train start context model from confident gene predictions.
/// Positive: seq_start of confident genes (high score, ATG, longest, long).
/// Negative: other in-frame ATG/GTG/TTG in same stop groups that weren't chosen.
pub fn train_start_model(seq: &[u8], orfs: &[Gene]) -> Option<StartModel> {
    // Positive: confident starts
    let confident: Vec<&Gene> = orfs.iter()
        .filter(|o| o.score > 0.50 && o.length >= 400 && o.is_atg() && o.is_longest)
        .collect();

    if confident.len() < 100 { return None; }

    let positives: Vec<usize> = confident.iter().map(|o| o.seq_start).collect();

    // Negative: other starts in same stop groups (not selected)
    let conf_stops: std::collections::HashSet<(bool, u32)> = confident.iter()
        .map(|o| (o.is_plus, o.stop_group)).collect();
    let conf_starts: std::collections::HashSet<usize> = confident.iter()
        .map(|o| o.seq_start).collect();

    let negatives: Vec<usize> = orfs.iter()
        .filter(|o| {
            conf_stops.contains(&(o.is_plus, o.stop_group))
            && !conf_starts.contains(&o.seq_start)
            && o.length >= 90
        })
        .map(|o| o.seq_start)
        .collect();

    if negatives.len() < 50 { return None; }

    // Train: 30 epochs, lr=0.01, L2=0.001
    StartModel::train(seq, &positives, &negatives, 30, 0.01, 0.001)
}

/// Score all ORFs with the trained start model.
pub fn score_start_model(seq: &[u8], orfs: &mut [Gene], model: &StartModel) {
    for orf in orfs.iter_mut() {
        orf.start_nn = model.score(seq, orf.seq_start);
    }
}

// ═══════════════════════════════════════════════════════════════
// N-terminal protein composition model
// ═══════════════════════════════════════════════════════════════

const NTERM_LEN: usize = 10; // score first 10 amino acids
const N_AA: usize = 20;

/// Translate a codon to amino acid index (0-19), or 20 for stop/unknown.
fn codon_to_aa(codon: &[u8]) -> usize {
    match codon {
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 0,  // A
        b"TGT" | b"TGC" => 1,                      // C
        b"GAT" | b"GAC" => 2,                      // D
        b"GAA" | b"GAG" => 3,                      // E
        b"TTT" | b"TTC" => 4,                      // F
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 5,  // G
        b"CAT" | b"CAC" => 6,                      // H
        b"ATT" | b"ATC" | b"ATA" => 7,            // I
        b"AAA" | b"AAG" => 8,                      // K
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 9, // L
        b"ATG" => 10,                               // M
        b"AAT" | b"AAC" => 11,                     // N
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 12,  // P
        b"CAA" | b"CAG" => 13,                     // Q
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 14, // R
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 15, // S
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 16,  // T
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 17,  // V
        b"TGG" => 18,                               // W
        b"TAT" | b"TAC" => 19,                     // Y
        _ => 20,                                     // stop or unknown
    }
}

/// N-terminal composition model: log-likelihood ratio at each position.
pub struct NtermModel {
    /// pwm[position][amino_acid] = log(freq_real / freq_background)
    pwm: Vec<[f64; N_AA]>,
}

impl NtermModel {
    /// Score the N-terminal amino acid sequence of an ORF.
    pub fn score(&self, seq: &[u8], start_pos: usize) -> f64 {
        let mut total = 0.0f64;
        let mut count = 0;
        for pos in 0..self.pwm.len() {
            let codon_start = start_pos + pos * 3;
            if codon_start + 3 > seq.len() { break; }
            let aa = codon_to_aa(&seq[codon_start..codon_start + 3]);
            if aa < N_AA {
                total += self.pwm[pos][aa];
                count += 1;
            }
        }
        if count > 0 { total / count as f64 } else { 0.0 }
    }
}

/// Train N-terminal model from confident gene starts.
/// Positive: N-terminal amino acids of confident genes.
/// Background: amino acid frequencies from random coding positions.
pub fn train_nterm_model(seq: &[u8], orfs: &[Gene]) -> Option<NtermModel> {
    let confident: Vec<&Gene> = orfs.iter()
        .filter(|o| o.score > 0.50 && o.length >= 400 && o.is_atg() && o.is_longest)
        .collect();
    if confident.len() < 100 { return None; }

    // Count amino acids at each N-terminal position
    let mut counts = vec![[0u32; N_AA]; NTERM_LEN];
    let mut n_seqs = 0u32;

    for o in &confident {
        let ss = o.seq_start;
        if ss + NTERM_LEN * 3 > seq.len() { continue; }
        let mut valid = true;
        for pos in 0..NTERM_LEN {
            let aa = codon_to_aa(&seq[ss + pos * 3..ss + pos * 3 + 3]);
            if aa >= N_AA { valid = false; break; }
            counts[pos][aa] += 1;
        }
        if valid { n_seqs += 1; }
    }
    if n_seqs < 80 { return None; }

    // Background: overall amino acid frequency from coding regions
    let mut bg = [0u64; N_AA];
    let mut bg_total = 0u64;
    for o in &confident {
        let ss = o.seq_start;
        let se = o.seq_end;
        let mut i = ss;
        while i + 3 <= se.min(seq.len()) {
            let aa = codon_to_aa(&seq[i..i + 3]);
            if aa < N_AA { bg[aa] += 1; bg_total += 1; }
            i += 3;
        }
    }
    if bg_total == 0 { return None; }

    // Build PWM: log(freq_position / freq_background)
    let mut pwm = vec![[0.0f64; N_AA]; NTERM_LEN];
    for pos in 0..NTERM_LEN {
        for aa in 0..N_AA {
            let freq = (counts[pos][aa] as f64 + 0.5) / (n_seqs as f64 + 0.5 * N_AA as f64);
            let bg_freq = (bg[aa] as f64 + 0.5) / (bg_total as f64 + 0.5 * N_AA as f64);
            pwm[pos][aa] = (freq / bg_freq).ln();
        }
    }

    Some(NtermModel { pwm })
}

/// Score all ORFs with the N-terminal model. Stores in stop_nn (reused as nterm_score).
pub fn score_nterm_model(seq: &[u8], orfs: &mut [Gene], model: &NtermModel) {
    for orf in orfs.iter_mut() {
        orf.stop_nn = model.score(seq, orf.seq_start);
    }
}

// ═══════════════════════════════════════════════════════════════
// Leader peptide detector (domain-based, not sequence-based)
// ═══════════════════════════════════════════════════════════════

/// Detect Rho-independent transcription terminator.
/// Structure: inverted repeat (stem-loop, stem ≥5bp) followed by ≥4 T's.
/// Scans a window starting at `pos` for `window` bp downstream.
pub fn detect_rho_independent_terminator(seq: &[u8], pos: usize, window: usize) -> bool {
    let end = (pos + window).min(seq.len());
    if end <= pos + 15 { return false; }
    let region = &seq[pos..end];

    // Look for inverted repeat (palindrome that forms stem-loop)
    // Scan for stems of length 5-15bp with loop of 3-8bp
    for stem_len in (5..=12).rev() {
        for start in 0..region.len().saturating_sub(stem_len * 2 + 3) {
            let left = &region[start..start + stem_len];
            // Check for complement in reverse after a loop
            for loop_len in 3..=8 {
                let right_start = start + stem_len + loop_len;
                if right_start + stem_len > region.len() { break; }
                let right = &region[right_start..right_start + stem_len];

                // Count complementary base pairs and GC content in stem
                let mut pairs = 0;
                let mut gc_pairs = 0;
                for k in 0..stem_len {
                    let l = left[k];
                    let r = right[stem_len - 1 - k];
                    if (l == b'G' && r == b'C') || (l == b'C' && r == b'G') {
                        pairs += 1;
                        gc_pairs += 1;
                    } else if (l == b'A' && r == b'T') || (l == b'T' && r == b'A') {
                        pairs += 1;
                    }
                }

                // Require ≥75% complementarity AND ≥40% GC pairs (real terminators are GC-rich)
                if pairs < (stem_len * 3 / 4).max(5) { continue; }
                if gc_pairs < (stem_len * 2 / 5).max(3) { continue; }

                // Check for polyT run after the stem-loop
                let poly_t_start = right_start + stem_len;
                let mut t_count = 0;
                for p in poly_t_start..(poly_t_start + 8).min(region.len()) {
                    if region[p] == b'T' { t_count += 1; }
                }

                // Require ≥5 T's in 8bp window
                if t_count >= 5 {
                    return true;
                }
            }
        }
    }
    false
}

/// Detect leader peptides at operon starts.
///
/// Strategy (operon-aware):
/// 1. Find operon starts: predicted genes with no same-strand upstream neighbor within 200bp
/// 2. Scan the region UPSTREAM of each operon-start gene (200-500bp window)
/// 3. Look for: short ORF + AA enrichment (domain) + attenuator (terminator)
///
/// Domain: abnormal enrichment of one amino acid (>20%, 4x expected).
/// thrL=Thr-rich, hisL=His-rich, pheL=Phe-rich — no sequence similarity.
pub fn find_leader_peptides(
    seq: &[u8],
    is_plus: bool,
    genome_len: usize,
    predicted_genes: &[Gene],
) -> Vec<Gene> {
    let mut leaders = Vec::new();
    let n = seq.len();

    // Find operon-start genes: first gene on this strand with no upstream same-strand neighbor
    let mut strand_genes: Vec<&Gene> = predicted_genes.iter()
        .filter(|g| g.is_plus == is_plus)
        .collect();
    strand_genes.sort_by_key(|g| g.start);

    let mut operon_starts: Vec<usize> = Vec::new(); // seq_start positions
    for (idx, g) in strand_genes.iter().enumerate() {
        let is_operon_start = if idx == 0 {
            true
        } else {
            let prev = strand_genes[idx - 1];
            // >200bp gap to previous same-strand gene = new operon
            if is_plus {
                g.start > prev.end + 200
            } else {
                prev.start > g.end + 200
            }
        };
        if is_operon_start {
            operon_starts.push(g.seq_start);
        }
    }

    // For each operon start, scan upstream region for leader peptide
    for &gene_seq_start in &operon_starts {
        // Scan 50-500bp upstream of the operon's first gene
        let scan_start = gene_seq_start.saturating_sub(500);
        let scan_end = gene_seq_start.saturating_sub(30); // at least 30bp gap to gene
        if scan_end <= scan_start { continue; }

        // Find all short ORFs in this upstream region
        for frame in 0..3usize {
            let mut i = scan_start + (frame - scan_start % 3 + 3) % 3; // align to frame
            while i + 3 <= scan_end {
                let codon = &seq[i..i+3];
                if codon != b"ATG" && codon != b"GTG" && codon != b"TTG" {
                    i += 3;
                    continue;
                }

                let orf_start = i;
                let mut aa_counts = [0u32; 20];
                let mut aa_total = 0u32;
                let mut j = i;
                let mut hit_stop = false;

                while j + 3 <= n {
                    let c = &seq[j..j+3];
                    if super::orf::is_stop(c) {
                        hit_stop = true;
                        break;
                    }
                    let aa = codon_to_aa(c);
                    if aa < 20 {
                        aa_counts[aa] += 1;
                        aa_total += 1;
                    }
                    j += 3;
                    if aa_total > 55 { break; }
                }

                if !hit_stop || aa_total < 10 || aa_total > 55 {
                    i += 3;
                    continue;
                }

                let orf_end = j + 3;
                let orf_len = orf_end - orf_start;

                // Must end before the downstream gene starts
                if orf_end > gene_seq_start {
                    i += 3;
                    continue;
                }

                // AA enrichment domain
                let max_aa = *aa_counts.iter().max().unwrap_or(&0);
                let enrichment = if aa_total > 1 {
                    max_aa as f64 / (aa_total - 1) as f64
                } else { 0.0 };

                // Also check for RUNS of enriched AA (≥3 consecutive)
                // Real leaders have: HHHHHH, TTTITITT, FFFAFFF — not scattered
                let enriched_aa_idx = aa_counts.iter().enumerate()
                    .max_by_key(|(_, &c)| c).map(|(i, _)| i).unwrap_or(20);
                let mut max_run = 0u32;
                let mut cur_run = 0u32;
                {
                    let mut p = orf_start;
                    while p + 3 <= j {
                        let aa = codon_to_aa(&seq[p..p+3]);
                        if aa == enriched_aa_idx { cur_run += 1; max_run = max_run.max(cur_run); }
                        else { cur_run = 0; }
                        p += 3;
                    }
                }

                if enrichment < 0.20 || max_aa < 3 || max_run < 2 {
                    i += 3;
                    continue;
                }

                // Attenuator (terminator) between ORF and gene
                let has_term = detect_rho_independent_terminator(seq, orf_end, 120);
                if !has_term {
                    i += 3;
                    continue;
                }

                // Not overlapping any predicted gene
                let (g_start, g_end) = if is_plus {
                    (orf_start + 1, orf_end)
                } else {
                    (genome_len - orf_end + 1, genome_len - orf_start)
                };
                let overlaps = predicted_genes.iter().any(|g| g.start < g_end && g.end > g_start);
                if overlaps {
                    i += 3;
                    continue;
                }

                let mut g = Gene::new();
                g.start = g_start.min(g_end);
                g.end = g_start.max(g_end);
                g.is_plus = is_plus;
                g.seq_start = orf_start;
                g.seq_end = orf_end;
                g.length = orf_len;
                g.start_codon = [seq[orf_start], seq[orf_start+1], seq[orf_start+2]];
                g.hex_avg = enrichment;
                g.score = 0.5 + enrichment;
                leaders.push(g);

                i += 3;
            }
        }
    }

    leaders
}

// ═══════════════════════════════════════════════════════════════
// Stop codon context model
// ═══════════════════════════════════════════════════════════════

/// Window for stop context: 30bp coding (upstream) + stop (3bp) + 30bp noncoding (downstream) = 63bp
const STOP_UP: usize = 30;
const STOP_DOWN: usize = 30;
const STOP_WIN: usize = STOP_UP + 3 + STOP_DOWN; // 63
const STOP_FEAT: usize = STOP_WIN * 4; // 252

/// Logistic regression model for stop codon context.
pub struct StopModel {
    weights: Vec<f64>,
    bias: f64,
}

impl StopModel {
    fn extract(seq: &[u8], stop_pos: usize) -> Option<Vec<f64>> {
        if stop_pos < STOP_UP { return None; }
        let end = stop_pos + 3 + STOP_DOWN;
        if end > seq.len() { return None; }

        let window = &seq[stop_pos - STOP_UP..end];
        let mut feat = vec![0.0f64; STOP_FEAT];
        for (i, &b) in window.iter().enumerate() {
            let n = nt4(b);
            if n > 3 { return None; }
            feat[i * 4 + n] = 1.0;
        }
        Some(feat)
    }

    fn sigmoid(x: f64) -> f64 {
        1.0 / (1.0 + (-x).exp())
    }

    pub fn score(&self, seq: &[u8], stop_pos: usize) -> f64 {
        match Self::extract(seq, stop_pos) {
            Some(feat) => {
                let z: f64 = feat.iter().zip(self.weights.iter())
                    .map(|(f, w)| f * w).sum::<f64>() + self.bias;
                Self::sigmoid(z)
            }
            None => 0.5,
        }
    }

    pub fn train(
        seq: &[u8],
        positives: &[usize],
        negatives: &[usize],
        n_epochs: usize,
        lr: f64,
        l2: f64,
    ) -> Option<StopModel> {
        if positives.len() < 50 || negatives.len() < 50 { return None; }

        let mut pos_feats: Vec<Vec<f64>> = Vec::new();
        let mut neg_feats: Vec<Vec<f64>> = Vec::new();
        for &p in positives {
            if let Some(f) = Self::extract(seq, p) { pos_feats.push(f); }
        }
        for &p in negatives {
            if let Some(f) = Self::extract(seq, p) { neg_feats.push(f); }
        }
        if pos_feats.len() < 40 || neg_feats.len() < 40 { return None; }

        let mut weights = vec![0.0f64; STOP_FEAT];
        let mut bias = 0.0f64;
        let n_per_class = pos_feats.len().min(neg_feats.len());

        for _epoch in 0..n_epochs {
            for feat in pos_feats.iter().take(n_per_class) {
                let z: f64 = feat.iter().zip(weights.iter()).map(|(f, w)| f * w).sum::<f64>() + bias;
                let err = Self::sigmoid(z) - 1.0;
                for (j, &f) in feat.iter().enumerate() {
                    weights[j] -= lr * (err * f + l2 * weights[j]);
                }
                bias -= lr * err;
            }
            for feat in neg_feats.iter().take(n_per_class) {
                let z: f64 = feat.iter().zip(weights.iter()).map(|(f, w)| f * w).sum::<f64>() + bias;
                let err = Self::sigmoid(z);
                for (j, &f) in feat.iter().enumerate() {
                    weights[j] -= lr * (err * f + l2 * weights[j]);
                }
                bias -= lr * err;
            }
        }
        Some(StopModel { weights, bias })
    }
}

/// Train stop context model.
/// Positive: stop codons of confident genes (real gene ends).
/// Negative: in-frame stop codons in intergenic regions (not gene ends).
pub fn train_stop_model(seq: &[u8], orfs: &[Gene]) -> Option<StopModel> {
    let confident: Vec<&Gene> = orfs.iter()
        .filter(|o| o.score > 0.50 && o.length >= 400 && o.is_longest)
        .collect();
    if confident.len() < 100 { return None; }

    // Positive: stop codon positions of confident genes
    // Stop codon is at seq_end - 3 (seq_end is exclusive)
    let positives: Vec<usize> = confident.iter()
        .filter(|o| o.seq_end >= 3)
        .map(|o| o.seq_end - 3)
        .collect();

    // Negative: random in-frame stop codons NOT at gene ends
    let gene_stops: std::collections::HashSet<usize> = confident.iter()
        .map(|o| o.seq_end - 3).collect();

    let mut negatives = Vec::new();
    // Scan for stop codons in frame 0 that aren't gene ends
    let mut i = 0;
    while i + 3 <= seq.len() {
        let codon = &seq[i..i+3];
        if super::orf::is_stop(codon) && !gene_stops.contains(&i) {
            negatives.push(i);
        }
        i += 3;
    }

    if negatives.len() < 50 { return None; }

    // Subsample negatives (there are far more non-gene stops than gene stops)
    if negatives.len() > positives.len() * 3 {
        // Take evenly spaced subset
        let step = negatives.len() / (positives.len() * 3);
        negatives = negatives.into_iter().step_by(step.max(1)).collect();
    }

    StopModel::train(seq, &positives, &negatives, 30, 0.01, 0.001)
}

/// Score all ORFs with the trained stop model.
pub fn score_stop_model(seq: &[u8], orfs: &mut [Gene], model: &StopModel) {
    for orf in orfs.iter_mut() {
        if orf.seq_end >= 3 {
            orf.stop_nn = model.score(seq, orf.seq_end - 3);
        }
    }
}
