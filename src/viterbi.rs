//! HMM gene finder with Viterbi decoding + self-training + RBS.
//!
//! Multi-pass approach:
//!   Pass 1: Basic hexamer Viterbi → get initial gene set
//!   Pass 2-3: Retrain hexamer on confident genes + intergenic model → refined
//!   Final: RBS-based start refinement (post-processing)

use crate::types::{Gene, HexModel, N_HEX};
use crate::io::hex_enc;
use crate::coding::{train_hex_from_set, merge_hex, blend_hex};
use crate::rbs::score_rbs_at;

/// Full HMM gene finder: self-training + intergenic model + RBS refinement.
pub fn hmm_gene_finder(
    genome: &[u8],
    rc: &[u8],
    hex_fwd: &HexModel,
    hex_rev: &HexModel,
    _intergenic_hex: Option<&HexModel>,
) -> Vec<Gene> {
    let n = genome.len();
    if n < 100 { return vec![]; }

    // === PASS 1: Basic Viterbi with initial hexamer ===
    let (pass1, _) = run_viterbi_both_strands(genome, rc, hex_fwd, hex_rev, n, None);

    // === Self-training: retrain on confident genes ===
    let confident: Vec<Gene> = pass1.iter()
        .filter(|g| g.length >= 400 && g.hex_avg > 0.05)
        .cloned().collect();

    if confident.len() < 100 {
        return pass1;
    }

    // Build intergenic model from pass 1 gaps
    let ig1 = build_intergenic_hex(genome, &pass1);

    // Self-training iteration 1: blend 0.25
    let conf_plus: Vec<Gene> = confident.iter().filter(|g| g.is_plus).cloned().collect();
    let conf_minus: Vec<Gene> = confident.iter().filter(|g| !g.is_plus).cloned().collect();
    let t2p = train_hex_from_set(genome, &conf_plus, 200);
    let t2m = train_hex_from_set(rc, &conf_minus, 200);
    let hex1 = match merge_hex(&t2p, &t2m) {
        Some(new_hex) => blend_hex(hex_fwd, &new_hex, 0.25),
        None => return pass1,
    };

    let (pass2, _) = run_viterbi_both_strands(genome, rc, &hex1, &hex1, n, Some(&ig1));

    // Self-training iteration 2: blend 0.20 on top of iteration 1
    let conf2: Vec<Gene> = pass2.iter()
        .filter(|g| g.length >= 400 && g.hex_avg > 0.05)
        .cloned().collect();

    let (mut best, mut all_cands) = if conf2.len() >= 100 {
        let c2p: Vec<Gene> = conf2.iter().filter(|g| g.is_plus).cloned().collect();
        let c2m: Vec<Gene> = conf2.iter().filter(|g| !g.is_plus).cloned().collect();
        let t3p = train_hex_from_set(genome, &c2p, 200);
        let t3m = train_hex_from_set(rc, &c2m, 200);
        match merge_hex(&t3p, &t3m) {
            Some(new_hex) => {
                let hex2 = blend_hex(&hex1, &new_hex, 0.20);
                let ig2 = build_intergenic_hex(genome, &pass2);
                let (p3, cands) = run_viterbi_both_strands(genome, rc, &hex2, &hex2, n, Some(&ig2));
                if p3.len() >= pass2.len() * 85 / 100 {
                    (p3, cands)
                } else {
                    // Re-run pass2 to get candidates
                    run_viterbi_both_strands(genome, rc, &hex1, &hex1, n, Some(&ig1))
                }
            }
            None => run_viterbi_both_strands(genome, rc, &hex1, &hex1, n, Some(&ig1)),
        }
    } else {
        run_viterbi_both_strands(genome, rc, &hex1, &hex1, n, Some(&ig1))
    };

    // Gap-fill: rescue high-quality candidates that didn't make DP selection
    let sel_set: std::collections::HashSet<(usize, usize, bool)> = best.iter()
        .map(|g| (g.start, g.end, g.is_plus)).collect();
    let non_selected: Vec<usize> = (0..all_cands.len())
        .filter(|&i| !sel_set.contains(&(all_cands[i].start, all_cands[i].end, all_cands[i].is_plus)))
        .collect();
    let mut sel_idx: Vec<usize> = (0..all_cands.len())
        .filter(|&i| sel_set.contains(&(all_cands[i].start, all_cands[i].end, all_cands[i].is_plus)))
        .collect();
    super::selection::gap_fill(&mut sel_idx, &all_cands, &non_selected, -0.1, 150);

    // Targeted gap rescue: find large gaps and insert good candidates
    super::selection::gap_targeted_rescue(&mut sel_idx, &all_cands, n, 500);

    best = sel_idx.iter().map(|&i| all_cands[i].clone()).collect();

    // RBS-based start refinement (post-processing)
    refine_hmm_starts(&mut best, genome, rc, n);

    best.sort_by_key(|g| g.start);
    best
}

/// Build intergenic hexamer model from gaps between genes.
fn build_intergenic_hex(genome: &[u8], genes: &[Gene]) -> HexModel {
    let mut counts = [0u32; N_HEX];
    let mut total = 0u64;
    let n = genome.len();

    // Mark positions covered by genes
    let mut covered = vec![false; n];
    for g in genes {
        let s = g.start.saturating_sub(1);
        let e = g.end.min(n);
        for p in s..e { covered[p] = true; }
    }

    // Count hexamers in uncovered regions
    for i in 0..n.saturating_sub(5) {
        if covered[i] { continue; }
        if let Some(idx) = hex_enc(&genome[i..i+6]) {
            counts[idx] += 1;
            total += 1;
        }
    }

    if total < 1000 { return [0.0; N_HEX]; }

    // Log-likelihood vs uniform
    let uniform = 1.0 / N_HEX as f64;
    let mut model = [0.0f64; N_HEX];
    for i in 0..N_HEX {
        let freq = (counts[i] as f64 + 0.1) / (total as f64 + 0.1 * N_HEX as f64);
        model[i] = (freq / uniform).ln();
    }
    model
}

/// Run Viterbi on both strands, all 3 frame offsets, merge results.
/// Returns (selected_genes, all_candidates).
fn run_viterbi_both_strands(
    genome: &[u8],
    rc: &[u8],
    hex_fwd: &HexModel,
    hex_rev: &HexModel,
    genome_len: usize,
    ig_hex: Option<&HexModel>,
) -> (Vec<Gene>, Vec<Gene>) {
    let mut all_genes = Vec::new();

    for offset in 0..3 {
        let mut fwd = viterbi_single_strand(genome, hex_fwd, offset, true, genome_len, ig_hex);
        let mut rev = viterbi_single_strand(rc, hex_rev, offset, false, genome_len, ig_hex);
        all_genes.append(&mut fwd);
        all_genes.append(&mut rev);
    }

    // Filter: require positive coding evidence (length-dependent)
    all_genes.retain(|g| {
        if g.length >= 600 { g.hex_avg > -0.4 }
        else if g.length >= 300 { g.hex_avg > -0.05 }
        else if g.length >= 150 { g.hex_avg > 0.10 }
        else { g.hex_avg > 0.20 }
    });

    // Overlap resolution via DP
    all_genes.sort_by(|a, b| {
        let sa = a.hex_avg * (a.length as f64).ln();
        let sb = b.hex_avg * (b.length as f64).ln();
        sb.partial_cmp(&sa).unwrap_or(std::cmp::Ordering::Equal)
    });

    for g in all_genes.iter_mut() {
        g.weight = g.hex_avg * (g.length as f64).ln().max(1.0);
        // Normalized score for rescue compatibility
        let hex_norm = ((g.hex_avg + 1.5) / 4.5).clamp(0.0, 1.0);
        let len_norm = 1.0 - (-(g.length as f64) / 300.0).exp();
        g.score = 0.4 * hex_norm + 0.3 * len_norm + 0.15 * g.start_type() + 0.15 * (g.rbs * 0.8);
        g.is_longest = true;
    }

    // Iterative DP with operon proximity boosting
    let sel_idx = super::selection::iterative_dp(&mut all_genes, 60, 3);
    let mut selected: Vec<Gene> = sel_idx.iter().map(|&i| all_genes[i].clone()).collect();
    selected.sort_by_key(|g| g.start);
    (selected, all_genes)
}

/// Viterbi for one strand, one frame offset.
fn viterbi_single_strand(
    seq: &[u8],
    hex: &HexModel,
    offset: usize,
    is_plus: bool,
    genome_len: usize,
    ig_hex: Option<&HexModel>,
) -> Vec<Gene> {
    let n = seq.len();
    let mut genes = Vec::new();

    // Precompute in-frame hexamer emission (coding)
    let mut emit: Vec<f64> = Vec::new();
    let mut positions: Vec<usize> = Vec::new();
    let mut ig_emit: Vec<f64> = Vec::new();
    let mut fb_emit: Vec<f64> = Vec::new(); // frame bias

    let mut i = offset;
    while i + 6 <= n {
        let score = if let Some(idx) = hex_enc(&seq[i..i+6]) {
            hex[idx]
        } else {
            -1.0
        };
        emit.push(score);
        positions.push(i);

        // Intergenic emission: average over 3 shifted hexamers to smooth
        let ig_score = if let Some(ig) = ig_hex {
            let mut s = 0.0f64;
            let mut c = 0u32;
            for shift in 0..3 {
                let p = i + shift;
                if p + 6 <= n {
                    if let Some(idx) = hex_enc(&seq[p..p+6]) {
                        s += ig[idx];
                        c += 1;
                    }
                }
            }
            if c > 0 { s / c as f64 } else { 0.0 }
        } else {
            0.0
        };
        ig_emit.push(ig_score);

        // Frame bias: in-frame hex vs average of other two frames
        let mut other_sum = 0.0f64;
        let mut other_n = 0u32;
        for off in 1..=2 {
            let p = if i >= off { i - off } else { i + 3 - off };
            if p + 6 <= n {
                if let Some(idx) = hex_enc(&seq[p..p+6]) {
                    other_sum += hex[idx];
                    other_n += 1;
                }
            }
        }
        let other_avg = if other_n > 0 { other_sum / other_n as f64 } else { 0.0 };
        fb_emit.push((score - other_avg).max(-1.0).min(2.0));

        i += 3;
    }

    if emit.is_empty() { return genes; }
    let m = emit.len();

    // DP arrays
    let mut dp_nc = vec![0.0f64; m + 1];
    let mut dp_cd = vec![f64::NEG_INFINITY; m + 1];
    let mut tr_nc = vec![0u8; m + 1];
    let mut tr_cd = vec![1u8; m + 1];

    // Use contrast scoring when intergenic model available
    let has_ig = ig_hex.is_some();

    for ci in 0..m {
        let pos = positions[ci];
        let codon = if pos + 3 <= n { &seq[pos..pos+3] } else { continue };
        let e = emit[ci];
        let ig_e = ig_emit[ci];

        let is_start = codon == b"ATG" || codon == b"GTG" || codon == b"TTG";
        let is_stop = codon == b"TAA" || codon == b"TAG" || codon == b"TGA";

        let fb = fb_emit[ci];
        // Coding emission: hex + intergenic contrast + frame bias
        let cd_emit = if has_ig {
            e + (e - ig_e).max(-0.5).min(1.0) * 0.15 + fb.max(0.0) * 0.12
        } else {
            e + fb.max(0.0) * 0.08
        };

        // === Coding state ===
        let stay_cd = if !is_stop {
            dp_cd[ci] + cd_emit
        } else {
            f64::NEG_INFINITY
        };

        let enter_cd = if is_start {
            let start_qual = if codon == b"ATG" { 1.0 } else if codon == b"GTG" { 0.5 } else { 0.3 };
            let rbs_bonus = score_rbs_at(seq, pos) * 0.8;
            dp_nc[ci] + cd_emit + start_qual + rbs_bonus
        } else {
            f64::NEG_INFINITY
        };

        if stay_cd >= enter_cd {
            dp_cd[ci + 1] = stay_cd;
            tr_cd[ci + 1] = 1;
        } else {
            dp_cd[ci + 1] = enter_cd;
            tr_cd[ci + 1] = 0;
        }

        // === Noncoding state ===
        // Noncoding emission: intergenic signal + anti-correlation with coding
        let nc_emit = if has_ig {
            0.02 + ig_e.max(-0.5).min(1.0) * 0.25 - e.max(0.0) * 0.12
        } else {
            0.1
        };

        let stay_nc = dp_nc[ci] + nc_emit;
        let exit_cd = if is_stop {
            dp_cd[ci] + 2.5
        } else {
            f64::NEG_INFINITY
        };

        if stay_nc >= exit_cd {
            dp_nc[ci + 1] = stay_nc;
            tr_nc[ci + 1] = 0;
        } else {
            dp_nc[ci + 1] = exit_cd;
            tr_nc[ci + 1] = 1;
        }
    }

    // Backtrack
    let mut state = if dp_nc[m] >= dp_cd[m] { 0u8 } else { 1 };
    let mut path = vec![0u8; m + 1];
    path[m] = state;
    for ci in (0..m).rev() {
        state = if state == 0 { tr_nc[ci + 1] } else { tr_cd[ci + 1] };
        path[ci] = state;
    }

    // Extract genes
    let mut in_gene = false;
    let mut gene_start_ci = 0;

    for ci in 0..=m {
        if !in_gene && ci < m && path[ci] == 1 {
            in_gene = true;
            gene_start_ci = ci;
        } else if in_gene && (ci == m || path[ci] == 0) {
            let seq_start = positions[gene_start_ci];
            let seq_end = if ci < m { positions[ci] } else { positions[m - 1] + 3 };
            let length = seq_end - seq_start;

            if length >= 90 {
                let (g_start, g_end) = if is_plus {
                    (seq_start + 1, seq_end)
                } else {
                    (genome_len - seq_end + 1, genome_len - seq_start)
                };

                let mut g = Gene::new();
                g.start = g_start.min(g_end);
                g.end = g_start.max(g_end);
                g.is_plus = is_plus;
                g.length = length;
                g.seq_start = seq_start;
                g.seq_end = seq_end;
                g.frame = offset as u8;
                if seq_start + 3 <= seq.len() {
                    g.start_codon = [seq[seq_start], seq[seq_start + 1], seq[seq_start + 2]];
                }
                let total_e: f64 = (gene_start_ci..ci.min(m)).map(|j| emit[j]).sum();
                let n_codons = (ci.min(m) - gene_start_ci).max(1);
                g.hex_avg = total_e / n_codons as f64;
                g.rbs = score_rbs_at(seq, seq_start);
                g.score = g.hex_avg;
                genes.push(g);
            }
            in_gene = false;
        }
    }

    genes
}

/// RBS-based start refinement for HMM-predicted genes.
/// For each gene, scan upstream in-frame for an alternative start
/// with better RBS + start_type + length consideration.
fn refine_hmm_starts(genes: &mut Vec<Gene>, genome: &[u8], rc: &[u8], genome_len: usize) {
    for g in genes.iter_mut() {
        let seq = if g.is_plus { genome } else { rc };
        let ss = g.seq_start;
        let se = g.seq_end;

        // Current start quality (includes length preference for current position)
        let cur_score = start_candidate_score(seq, ss, se);

        // Scan upstream: same frame, up to 180bp (60 codons)
        let mut best_pos = ss;
        let mut best_score = cur_score;
        let mut pos = ss;

        loop {
            if pos < 3 { break; }
            pos -= 3;
            if ss - pos > 180 { break; }
            if pos + 3 > seq.len() { break; }
            let codon = &seq[pos..pos+3];

            // Hit a stop codon → can't go further upstream
            if codon == b"TAA" || codon == b"TAG" || codon == b"TGA" {
                break;
            }

            if codon == b"ATG" || codon == b"GTG" || codon == b"TTG" {
                let score = start_candidate_score(seq, pos, se);
                // Require clear improvement to change start
                if score > best_score + 0.03 {
                    best_score = score;
                    best_pos = pos;
                }
            }
        }

        // Also scan downstream (truncation): find shorter start with much better RBS
        pos = ss;
        let max_downstream = (se - ss).min(120); // don't truncate more than 120bp
        let mut d_pos = ss + 3;
        while d_pos - ss < max_downstream && d_pos + 3 <= se {
            if d_pos + 3 > seq.len() { break; }
            let codon = &seq[d_pos..d_pos+3];
            if codon == b"TAA" || codon == b"TAG" || codon == b"TGA" { break; }
            if codon == b"ATG" || codon == b"GTG" || codon == b"TTG" {
                let score = start_candidate_score(seq, d_pos, se);
                // Higher threshold for truncation (losing codons must be justified)
                if score > best_score + 0.12 {
                    best_score = score;
                    best_pos = d_pos;
                }
            }
            d_pos += 3;
        }

        if best_pos != ss {
            let new_ss = best_pos;
            let new_length = se - new_ss;
            if new_length < 90 { continue; } // don't create too-short genes

            g.seq_start = new_ss;
            g.length = new_length;
            if new_ss + 3 <= seq.len() {
                g.start_codon = [seq[new_ss], seq[new_ss + 1], seq[new_ss + 2]];
            }

            if g.is_plus {
                g.start = new_ss + 1;
            } else {
                g.end = genome_len - new_ss;
            }

            g.rbs = score_rbs_at(seq, new_ss);
        }
    }
}

/// Score a start codon candidate by RBS + start_type + upstream AT content.
fn start_candidate_score(seq: &[u8], pos: usize, stop_pos: usize) -> f64 {
    if pos + 3 > seq.len() { return 0.0; }
    let codon = &seq[pos..pos+3];
    let start_type = if codon == b"ATG" { 1.0 }
        else if codon == b"GTG" { 0.55 }
        else if codon == b"TTG" { 0.35 }
        else { 0.0 };
    let rbs = score_rbs_at(seq, pos);

    // Upstream AT content (SD motifs are AT-rich context)
    let up_at = if pos >= 20 {
        let region = &seq[pos-20..pos];
        region.iter().filter(|&&c| c == b'A' || c == b'T').count() as f64 / 20.0
    } else {
        0.5
    };
    let at_bonus = ((up_at - 0.4) * 0.3).max(0.0).min(0.1);

    // Small length preference: longer genes are slightly preferred
    let gene_len = stop_pos - pos;
    let len_bonus = (1.0 - (-(gene_len as f64) / 500.0).exp()) * 0.08;

    rbs * 0.48 + start_type * 0.35 + at_bonus + len_bonus
}

/// Compute viterbi_frac for existing pipeline genes.
pub fn compute_viterbi_frac(
    seq: &[u8],
    _rc: &[u8],
    hex: &HexModel,
    genes: &mut [Gene],
) {
    let n = seq.len();
    let mut coding = vec![false; n];
    for offset in 0..3 {
        let mut cum = 0.0f64;
        let mut in_orf = false;
        let mut orf_start = 0;
        let mut i = offset;
        while i + 3 <= n {
            let codon = &seq[i..i+3];
            let is_start = codon == b"ATG" || codon == b"GTG" || codon == b"TTG";
            let is_stop = codon == b"TAA" || codon == b"TAG" || codon == b"TGA";
            if !in_orf && is_start {
                in_orf = true;
                orf_start = i;
                cum = 0.0;
            }
            if in_orf {
                if i + 6 <= n {
                    if let Some(idx) = hex_enc(&seq[i..i+6]) {
                        cum += hex[idx];
                    }
                }
                if is_stop {
                    if cum > 0.0 && i - orf_start >= 87 {
                        for p in orf_start..i+3 { coding[p] = true; }
                    }
                    in_orf = false;
                }
            }
            i += 3;
        }
    }
    for gene in genes.iter_mut() {
        let start = gene.start.saturating_sub(1);
        let end = gene.end.min(n);
        if start < end {
            let c = coding[start..end].iter().filter(|&&x| x).count();
            gene.viterbi_frac = c as f64 / (end - start) as f64;
        }
    }
}
