use crate::types::{Gene, HexModel};
use crate::io::hex_enc;

/// Probability of NOT hitting a stop codon per codon, given genome GC.
/// Prodigal formula (node.c:575): for translation table 11 (standard):
///   P(stop) = (1-gc)^2 * gc/4 + (1-gc)^3 / 8
///   no_stop = 1 - P(stop)
/// At GC=0.30: no_stop ≈ 0.974 (stops frequent, length meaningful)
/// At GC=0.50: no_stop ≈ 0.984
/// At GC=0.70: no_stop ≈ 0.993 (stops rare, length less meaningful)
fn no_stop_probability(gc: f64) -> f64 {
    let gc = gc.clamp(0.15, 0.90);
    let at = 1.0 - gc;
    let p_stop = at * at * gc / 4.0 + at * at * at / 8.0;
    1.0 - p_stop
}

/// Composite score.
pub fn composite_score(g: &Gene, gc3_target: f64) -> f64 {
    let len = g.length as f64;
    let len_norm = 1.0 - (-len / 300.0).exp();
    let st = g.start_type();
    let hex_norm = ((g.hex_avg + 1.5) / 4.5).clamp(0.0, 1.0);
    let fb_norm = (g.frame_bias / 3.0).clamp(0.0, 1.0);
    let cov = g.hex_cov;
    let mono_norm = ((g.mono + 0.5) / 2.0).clamp(0.0, 1.0);
    let gc3_dev = (g.gc3 - gc3_target).abs();
    let gc3_norm = (1.0 - gc3_dev * 4.0).max(0.0);
    let edge_norm = ((g.edge + 0.5) / 3.0).clamp(0.0, 1.0);

    // rbs_combined: always includes rbs_pwm floor (matches monolith exactly)
    let rbs_pwm_norm = ((g.rbs_pwm + 3.0) / 13.0).clamp(0.0, 1.0);
    let rbs_combined = g.rbs.max(rbs_pwm_norm * 0.85);

    // Upstream AT bonus: promoter regions are AT-rich relative to genome background.
    // Normalize by genome GC: compare upstream_at vs expected AT = (1 - gc3_target).
    let genome_at = 1.0 - gc3_target.clamp(0.2, 0.85);
    let at_bonus = ((g.upstream_at - genome_at) * 0.5).clamp(0.0, 0.10);
    let longest_bonus = if g.is_longest { 0.06 } else { 0.0 };

    // ML-informed weights (logistic regression on 192K ORFs):
    // Top signals: hex_avg (2.14), length (2.40), upstream_at (0.82), is_longest (0.74), rbs (0.63)
    // Weak/zero: mono (0), edge (0), shadow_pen (0)
    // Kept but reduced: frame_bias, gc3, cov
    let mut score = 0.28 * hex_norm        // ML: 2.14 — dominant
        + 0.04 * fb_norm                   // ML: 0.38 — much less important than we thought
        + 0.14 * rbs_combined              // ML: 0.63
        + 0.15 * len_norm                  // ML: 2.40 — was 0.08, massively underweighted!
        + 0.05 * st                        // start type
        + 0.03 * gc3_norm                  // ML: 0.53
        + 0.01 * edge_norm                 // ML: 0 — nearly useless
        + 0.04 * cov                       // ML: 0.24
        + 0.00 * mono_norm                 // ML: 0 — useless, removed
        + 0.06 * at_bonus                  // ML: 0.82 — was 0.02, hugely underweighted!
        + 0.08 * longest_bonus;            // ML: 0.74 — was 0.05

    // TIR composite — DISABLED: F1 0.939→0.937, linear weights already capture interaction
    // let tir_quality = rbs_combined * st * (0.5 + g.upstream_at);

    // RBS rescue for clearly coding ORFs
    if rbs_combined < 0.10 && hex_norm > 0.55 && g.is_atg() {
        score += 0.04;
    }
    // Leaderless gene bonus (cycle 20): promoter compensates for missing SD
    if rbs_combined < 0.20 && g.leaderless > 0.8 && g.is_atg() {
        score += 0.05;
    }
    // Start context bonus (cycle 25): Kozak-like context around ATG
    score += 0.02 * (g.start_ctx / 0.38).min(1.0);
    // GC3 bias bonus (cycle 15): wobble position bias indicates coding
    score += 0.02 * (g.gc3_bias * 8.0).min(1.0);
    // Start context neural model: used only in keep_best_starts for start selection,
    // NOT in composite scoring (adding it here adds FP without improving F1).

    score
}

/// GC3 bias: |GC3 - GC12| (cycle 15 idea, survived crash).
/// Wobble position (3rd codon pos) has different GC than positions 1+2.
/// Coding regions show larger GC3 bias than non-coding.
pub fn compute_gc3_bias(seq: &[u8], ss: usize, se: usize) -> f64 {
    let s = &seq[ss..se];
    let mut gc = [0u32; 3];
    let mut tot = [0u32; 3];
    let mut i = 0;
    while i + 2 < s.len() {
        for p in 0..3 {
            tot[p] += 1;
            if s[i + p] == b'G' || s[i + p] == b'C' { gc[p] += 1; }
        }
        i += 3;
    }
    if tot[0] + tot[1] == 0 || tot[2] == 0 { return 0.0; }
    let gc3 = gc[2] as f64 / tot[2] as f64;
    let gc12 = (gc[0] + gc[1]) as f64 / (tot[0] + tot[1]) as f64;
    (gc3 - gc12).abs()
}

pub fn compute_gc3(seq: &[u8], ss: usize, se: usize) -> f64 {
    let s = &seq[ss..se];
    let mut gc = 0u32;
    let mut total = 0u32;
    let mut i = 2;
    while i < s.len().saturating_sub(2) {
        match s[i] {
            b'G' | b'C' => { gc += 1; total += 1; }
            b'A' | b'T' => { total += 1; }
            _ => {}
        }
        i += 3;
    }
    if total > 0 { gc as f64 / total as f64 } else { 0.5 }
}

/// Upstream in-frame coding potential.
/// Scores the region upstream of start_pos IN THE SAME READING FRAME.
/// If this score is high → start_pos is likely an internal start (gene continues upstream).
/// If low → start_pos is likely the real start (upstream is non-coding in this frame).
/// Scans up to `max_upstream` bp, stops at in-frame stop codon.
pub fn upstream_inframe_coding(seq: &[u8], start_pos: usize, model: &HexModel, max_upstream: usize) -> f64 {
    let mut total = 0.0f64;
    let mut n = 0u32;

    // Go upstream in steps of 3 (same frame)
    let limit = start_pos.saturating_sub(max_upstream);
    let mut pos = start_pos;

    loop {
        if pos < 3 { break; }
        pos -= 3;
        if pos < limit { break; }

        // Stop at in-frame stop codon (this is the real ORF boundary)
        if pos + 3 <= seq.len() {
            let codon = &seq[pos..pos+3];
            if codon == b"TAA" || codon == b"TAG" || codon == b"TGA" {
                break;
            }
        }

        // Score hexamer at this position
        if pos + 6 <= seq.len() {
            if let Some(idx) = hex_enc(&seq[pos..pos+6]) {
                total += model[idx];
                n += 1;
            }
        }
    }

    if n > 0 { total / n as f64 } else { 0.0 }
}

/// Edge coding score matching monolith's score_coding_edge.
/// Upstream iterates every position (all frames), downstream every 3 (in-frame).
/// Does NOT clamp to max(0.0) — negative values are informative.
pub fn edge_coding_score(seq: &[u8], start_pos: usize, model: &HexModel, window: usize) -> f64 {
    let up_start = start_pos.saturating_sub(window);

    // Upstream: iterate every position (checking all reading frames)
    let mut up_total = 0.0f64;
    let mut up_n = 0u32;
    if start_pos >= 5 {
        for i in up_start..start_pos.saturating_sub(5) {
            if i + 6 <= seq.len() {
                if let Some(idx) = hex_enc(&seq[i..i+6]) {
                    up_total += model[idx];
                    up_n += 1;
                }
            }
        }
    }

    // Downstream: iterate every 3 (in-frame)
    let mut down_total = 0.0f64;
    let mut down_n = 0u32;
    let end_pos = (start_pos + window).min(seq.len().saturating_sub(5));
    let mut i = start_pos;
    while i < end_pos {
        if i + 6 <= seq.len() {
            if let Some(idx) = hex_enc(&seq[i..i+6]) {
                down_total += model[idx];
                down_n += 1;
            }
        }
        i += 3;
    }

    let up_avg = if up_n > 0 { up_total / up_n as f64 } else { 0.0 };
    let down_avg = if down_n > 0 { down_total / down_n as f64 } else { 0.0 };
    down_avg - up_avg
}
