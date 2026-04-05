//! GC Frame Plot: zero-knowledge coding signal for bootstrap.
//!
//! In coding regions, the 3rd codon position (wobble) has systematically
//! different GC content than positions 1 and 2. This provides a coding
//! signal without any trained model.
//!
//! Key insight: for an ORF starting at genome position X, the wobble
//! positions are at X+2, X+5, X+8, ... (i.e. genome frame (X+2)%3).
//! A real gene shows consistent GC bias at its wobble frame.
//! A random ORF does not — the bias is at an arbitrary frame.
//!
//! Based on Prodigal's approach (Hyatt et al 2010).

use crate::types::Gene;

const WINDOW: usize = 120;

/// For each position in the sequence, compute which of the 3 genome
/// frames (0,1,2) has the highest GC content in a surrounding window.
/// Returns array of winning frame indices.
pub fn calc_gc_frame_bias(seq: &[u8]) -> Vec<u8> {
    let slen = seq.len();
    if slen < WINDOW {
        return vec![0; slen];
    }

    let mut result = vec![0u8; slen];
    let half = WINDOW / 2;

    for center in 0..slen {
        let start = center.saturating_sub(half);
        let end = (center + half).min(slen);

        let mut gc = [0u32; 3];
        let mut total = [0u32; 3];

        for i in start..end {
            let frame = i % 3;
            total[frame] += 1;
            if seq[i] == b'G' || seq[i] == b'C' {
                gc[frame] += 1;
            }
        }

        let mut best_frame = 0u8;
        let mut best_frac = 0.0f64;
        for f in 0..3 {
            if total[f] > 0 {
                let frac = gc[f] as f64 / total[f] as f64;
                if frac > best_frac {
                    best_frac = frac;
                    best_frame = f as u8;
                }
            }
        }
        result[center] = best_frame;
    }

    result
}

/// Compute genome-wide bias: which frame tends to have highest GC.
/// Returns bias[3] weights normalized so they sum to 3.0.
/// Computed from all long ORFs (>= min_len).
pub fn compute_genome_bias(gc_bias: &[u8], orfs: &[Gene], min_len: usize) -> [f64; 3] {
    let mut frame_weight = [0.0f64; 3];

    for orf in orfs {
        if orf.length < min_len || !orf.is_longest {
            continue;
        }

        let mut counts = [0u32; 3];
        let mut total = 0u32;

        let mut pos = orf.seq_start;
        while pos < orf.seq_end.min(gc_bias.len()) {
            let winning = gc_bias[pos] as usize;
            if winning < 3 {
                counts[winning] += 1;
            }
            total += 1;
            pos += 1;
        }

        if total > 0 {
            let len_factor = orf.length as f64 / 1000.0;
            for j in 0..3 {
                frame_weight[j] += 3.0 * counts[j] as f64 / total as f64 * len_factor;
            }
        }
    }

    // Normalize to sum = 3.0
    let sum: f64 = frame_weight.iter().sum();
    if sum > 0.0 {
        for j in 0..3 {
            frame_weight[j] = frame_weight[j] * 3.0 / sum;
        }
    } else {
        frame_weight = [1.0, 1.0, 1.0];
    }

    frame_weight
}

/// Frame-aware GC score for a single ORF.
///
/// For an ORF starting at seq_start, the wobble position is at
/// (seq_start + 2) % 3 in genome coordinates. We check what fraction
/// of positions show the EXPECTED wobble frame as the GC-dominant frame.
///
/// Returns: weighted score using genome-wide bias weights.
/// High score = this ORF's GC frame pattern is consistent with coding.
pub fn score_gc_frame_aware(gc_bias: &[u8], orf: &Gene, genome_bias: &[f64; 3]) -> f64 {
    if orf.seq_end <= orf.seq_start || orf.seq_end > gc_bias.len() {
        return 0.0;
    }

    // Per-ORF: count how often each genome frame wins, weighted by position in codon
    let mut gc_score = [0.0f64; 3];
    let mut total = 0u32;

    let mut pos = orf.seq_start;
    while pos + 2 < orf.seq_end && pos + 2 < gc_bias.len() {
        // This codon spans genome positions pos, pos+1, pos+2
        // The winning GC frame at each position
        for offset in 0..3 {
            let p = pos + offset;
            if p < gc_bias.len() {
                let winning = gc_bias[p] as usize;
                if winning < 3 {
                    gc_score[winning] += 1.0;
                    total += 1;
                }
            }
        }
        pos += 3;
    }

    if total == 0 {
        return 0.0;
    }

    // Normalize to fractions (should sum to ~1.0)
    for j in 0..3 {
        gc_score[j] /= total as f64;
    }

    // Weighted score: dot product with genome-wide bias
    // Real genes: their gc_score distribution matches genome_bias → high dot product
    // Random ORFs: gc_score is at a random frame → lower dot product
    let score = genome_bias[0] * gc_score[0]
              + genome_bias[1] * gc_score[1]
              + genome_bias[2] * gc_score[2];

    // Normalize: uniform distribution gives score = 1.0 (bias·[1/3,1/3,1/3] = 1.0)
    // Perfect match gives score up to ~3.0 (if one frame has bias=3.0 and gc_score=1.0)
    // Map to 0..1 range: (score - 1.0) / 2.0
    ((score - 1.0) / 2.0).clamp(0.0, 1.0)
}

/// Score all ORFs on one strand using frame-aware GC bias.
/// Stores result in orf.adj_bonus.
pub fn score_all_gc_frame(seq: &[u8], orfs: &mut [Gene]) {
    let gc_bias = calc_gc_frame_bias(seq);
    let genome_bias = compute_genome_bias(&gc_bias, orfs, 600);
    for orf in orfs.iter_mut() {
        orf.adj_bonus = score_gc_frame_aware(&gc_bias, orf, &genome_bias);
    }
}

/// Initial DP using only length + GC frame score (no hexamers).
/// Selects a clean set of likely-real genes for hexamer training.
/// Returns indices into the orfs slice.
pub fn gc_frame_select(orfs: &[Gene], gc3_target: f64) -> Vec<usize> {
    // Score each ORF: length significance × GC frame consistency
    let gc = gc3_target.clamp(0.15, 0.90);
    let at = 1.0 - gc;
    let p_stop = at * at * gc / 4.0 + at * at * at / 8.0;
    let ns = 1.0 - p_stop;

    let mut scored: Vec<(usize, f64)> = Vec::new();

    for (i, orf) in orfs.iter().enumerate() {
        if !orf.is_longest {
            continue;
        }
        let n_codons = orf.length as f64 / 3.0;
        let len_sig = 1.0 - ns.powf(n_codons); // how unlikely this length is by chance
        let gc_score = orf.adj_bonus; // frame-aware GC score

        // Combined: both length AND GC frame must be significant
        let combined = len_sig * 0.5 + gc_score * 0.5;

        if combined > 0.25 && orf.length >= 300 {
            scored.push((i, combined));
        }
    }

    // Sort by score descending
    scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    // Simple greedy non-overlapping selection
    let mut selected: Vec<usize> = Vec::new();
    let mut used: Vec<(usize, usize)> = Vec::new(); // (start, end) of selected genes

    for (idx, _score) in &scored {
        let orf = &orfs[*idx];
        let overlaps = used.iter().any(|&(s, e)| {
            orf.start < e && orf.end > s
        });
        if !overlaps {
            selected.push(*idx);
            used.push((orf.start, orf.end));
        }
    }

    selected
}
