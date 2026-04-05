//! GC Frame Plot: zero-knowledge coding signal for bootstrap.
//!
//! In coding regions, the 3rd codon position (wobble) has different GC
//! content than positions 1 and 2. This asymmetry provides a coding signal
//! without any trained model:
//! - High-GC genomes: position 2 (3rd codon pos in frame 0) tends to have highest GC
//! - Low-GC genomes: position 0 (1st codon pos) tends to have highest GC
//!
//! Based on Prodigal's approach (Hyatt et al 2010).

use crate::types::Gene;

const WINDOW: usize = 120;

/// For each position in the sequence, compute which codon frame (0,1,2)
/// has the highest GC content in a surrounding window.
pub fn calc_gc_frame_bias(seq: &[u8]) -> Vec<u8> {
    let slen = seq.len();
    if slen < WINDOW {
        return vec![0; slen];
    }

    let mut result = vec![0u8; slen];

    // Sliding window: count GC per frame position
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

        // Which frame has highest GC fraction?
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

/// Score an ORF by GC frame consistency: what fraction of its codons
/// have the expected frame as highest-GC.
///
/// For a real gene in frame 0, the 3rd codon position (frame offset 2)
/// should consistently win. Returns a score 0..1 where higher = more
/// consistent with coding.
pub fn score_gc_frame(gc_bias: &[u8], orf: &Gene) -> f64 {
    if orf.seq_end <= orf.seq_start || orf.seq_end > gc_bias.len() {
        return 0.0;
    }

    let mut frame_counts = [0u32; 3];
    let mut total = 0u32;

    // Walk codon by codon through the ORF
    let mut pos = orf.seq_start;
    while pos + 2 < orf.seq_end && pos < gc_bias.len() {
        // For each codon position, check which frame "wins" at that spot
        for offset in 0..3 {
            if pos + offset < gc_bias.len() {
                let winning_frame = gc_bias[pos + offset] as usize;
                frame_counts[winning_frame] += 1;
                total += 1;
            }
        }
        pos += 3;
    }

    if total == 0 {
        return 0.0;
    }

    // Compute genome-wide bias: which frame tends to have highest GC overall
    // For real coding ORFs, one frame should dominate much more than 1/3
    let fracs: [f64; 3] = [
        frame_counts[0] as f64 / total as f64,
        frame_counts[1] as f64 / total as f64,
        frame_counts[2] as f64 / total as f64,
    ];

    // Score = how far the dominant frame deviates from uniform (1/3)
    // Max deviation: 1.0 (one frame always wins) → score = 1.0
    // Uniform: each ~0.33 → score = 0.0
    let max_frac = fracs.iter().cloned().fold(0.0f64, f64::max);
    let score = ((max_frac - 0.333) / 0.667).clamp(0.0, 1.0);

    score
}

/// Compute GC frame bias scores for all ORFs on one strand.
/// Returns a vec of scores (0..1) parallel to the input orfs.
pub fn score_all_gc_frame(seq: &[u8], orfs: &mut [Gene]) {
    let gc_bias = calc_gc_frame_bias(seq);
    for orf in orfs.iter_mut() {
        let gc_score = score_gc_frame(&gc_bias, orf);
        // Use gc_frame as a bonus to the initial score before hexamer training
        // Store in adj_bonus temporarily (it gets overwritten later anyway)
        orf.adj_bonus = gc_score;
    }
}
