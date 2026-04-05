//! Rho-independent transcription terminator detection.
//!
//! Algorithm cannibalized from TransTermHP (Kingsford et al. 2007):
//! - DP-based hairpin scoring with mismatches and gaps
//! - Poly-U tail scoring with exponential decay
//! - U-window prefilter for speed
//!
//! Energy model (Ermolaeva/TransTermHP):
//!   G-C pair: -2.3, A-U pair: -0.9, G-U wobble: +1.3
//!   Mismatch: +3.5, Gap (bulge): +6.0
//!   Loop penalty: linear (loop_size - 2)

/// A detected Rho-independent terminator.
#[derive(Clone)]
pub struct Terminator {
    pub start: usize,     // genomic start (1-based)
    pub end: usize,       // genomic end (1-based)
    pub strand: bool,     // true = plus
    pub hairpin_score: f64,
    pub tail_score: f64,
    pub confidence: f64,
}

// Energy constants (TransTermHP defaults)
const GC: f64 = -2.3;
const AU: f64 = -0.9;
const GU: f64 = 1.3;
const MM: f64 = 3.5;
const GAP: f64 = 6.0;

const MIN_STEM: usize = 4;
const MIN_LOOP: usize = 3;
const MAX_LOOP: usize = 13;
const MAX_HP: usize = 59;  // max hairpin extent

// Prefilter: require ≥3 T's in window of 6 (TransTermHP UWINDOW)
const UWINDOW_SIZE: usize = 6;
const UWINDOW_REQUIRE: usize = 3;

// Thresholds
const HP_CUTOFF: f64 = -2.0;
const TAIL_CUTOFF: f64 = -2.5;

/// Loop penalty: linear, indexed by loop size.
#[inline]
fn loop_penalty(len: usize) -> f64 {
    if len < MIN_LOOP { return 1000.0; }
    if len > MAX_LOOP { return 1000.0; }
    (len as f64 - 2.0).max(0.0)
}

/// Score a base pair on the forward strand.
#[inline]
fn pair_score(a: u8, b: u8) -> f64 {
    match (a, b) {
        (b'G', b'C') | (b'C', b'G') => GC,
        (b'A', b'T') | (b'T', b'A') => AU,
        (b'G', b'T') | (b'T', b'G') => GU,
        _ => MM,
    }
}

/// Tail score: exponential decay, T gets 0.9x, non-T gets 0.6x.
/// Scans 15bp after the hairpin end.
fn tail_score_at(seq: &[u8], pos: usize) -> f64 {
    let mut sum = 0.0f64;
    let mut prev = 1.0f64;
    let end = (pos + 15).min(seq.len());
    for i in pos..end {
        prev *= if seq[i] == b'T' { 0.9 } else { 0.6 };
        sum += prev;
    }
    -sum
}

/// DP hairpin finder for a window of the sequence.
/// Returns (best_hp_score, stem_start_offset, stem_end_offset) relative to window start.
fn dp_best_hairpin(seq: &[u8], window_start: usize, window_len: usize) -> Option<(f64, usize, usize)> {
    let n = window_len.min(MAX_HP);
    if n < MIN_STEM * 2 + MIN_LOOP { return None; }

    let w = &seq[window_start..window_start + n];

    // DP table: dp[i][j] = best hairpin score for subsequence w[i..=j]
    // Using flat array for speed
    let mut dp = vec![1000.0f64; n * n];
    let mut reason = vec![0u8; n * n]; // 0=loop, 1=match/mm, 2=igap, 3=jgap

    #[inline]
    fn idx(i: usize, j: usize, n: usize) -> usize { i * n + j }

    // Fill DP: small spans first (inner loops), then larger (extending stems outward)
    // dp[i][j] = best hairpin score for w[i..=j]
    // Base: loop_penalty for small spans
    // Extension: pair w[i] with w[j], add inner dp[i+1][j-1]
    for span in MIN_LOOP..n {
        for i in 0..n.saturating_sub(span) {
            let j = i + span;
            let k = idx(i, j, n);

            // Option 1: entire span is a loop
            dp[k] = loop_penalty(span + 1);
            reason[k] = 0;

            // Option 2: extend stem by pairing w[i] with w[j]
            if span >= MIN_LOOP + 2 {
                let pair_e = pair_score(w[i], w[j]);
                let inner = dp[idx(i + 1, j - 1, n)];
                let extend = pair_e + inner;
                if extend <= dp[k] {
                    dp[k] = extend;
                    reason[k] = 1;
                }
            }

            // Option 3: gap on left (bulge: skip w[i], pair w[i+1] with w[j])
            if span >= MIN_LOOP + 3 && i + 2 <= j - 1 {
                let gap_l = GAP + dp[idx(i + 2, j - 1, n)];
                if gap_l <= dp[k] {
                    dp[k] = gap_l;
                    reason[k] = 2;
                }
            }

            // Option 4: gap on right (bulge: pair w[i] with w[j-1], skip w[j])
            if span >= MIN_LOOP + 3 && i + 1 <= j - 2 {
                let gap_r = GAP + dp[idx(i + 1, j - 2, n)];
                if gap_r <= dp[k] {
                    dp[k] = gap_r;
                    reason[k] = 3;
                }
            }
        }
    }

    // Find best hairpin anchored at position 0
    let mut best_score = 0.0f64; // must be negative
    let mut best_j = 0;
    for j in (MIN_STEM * 2 + MIN_LOOP - 1)..n {
        let s = dp[idx(0, j, n)];
        if s < best_score {
            best_score = s;
            best_j = j;
        }
    }

    if best_score < HP_CUTOFF {
        Some((best_score, 0, best_j + 1))
    } else {
        None
    }
}

/// Find all Rho-independent terminators on one strand.
pub fn find_terminators(seq: &[u8], is_plus: bool, genome_len: usize) -> Vec<Terminator> {
    let mut terminators = Vec::new();
    let n = seq.len();
    if n < MAX_HP + 15 { return terminators; }

    // Precompute T-count sliding window
    let mut t_window = 0u32;
    for i in 0..UWINDOW_SIZE.min(n) {
        if seq[i] == b'T' { t_window += 1; }
    }

    let mut pos = 0usize;
    let end_pos = n.saturating_sub(MAX_HP + 15);

    while pos < end_pos {
        // TransTermHP prefilter: ≥3 T's in 6bp window after current position
        // AND current position is NOT T (the hairpin end is not itself in the tail)
        if t_window >= UWINDOW_REQUIRE as u32 && seq[pos] != b'T' {
            let wlen = MAX_HP.min(n - pos - 15);
            if let Some((hp_score, _hp_start, hp_end)) = dp_best_hairpin(seq, pos, wlen) {
                let ts = tail_score_at(seq, pos + hp_end);
                if ts < TAIL_CUTOFF {
                    let (g_start, g_end) = if is_plus {
                        (pos + 1, pos + hp_end + 15)
                    } else {
                        let s = genome_len - (pos + hp_end + 15) + 1;
                        let e = genome_len - pos;
                        (s.min(e), s.max(e))
                    };

                    let confidence = (-hp_score * 8.0 + (-ts) * 12.0).min(100.0);

                    terminators.push(Terminator {
                        start: g_start,
                        end: g_end,
                        strand: is_plus,
                        hairpin_score: hp_score,
                        tail_score: ts,
                        confidence,
                    });

                    // Skip past this terminator to avoid overlapping hits
                    let skip = hp_end.max(5);
                    for _ in 0..skip {
                        if pos + UWINDOW_SIZE < n && seq[pos] == b'T' { t_window -= 1; }
                        pos += 1;
                        if pos + UWINDOW_SIZE < n && seq[pos + UWINDOW_SIZE] == b'T' { t_window += 1; }
                    }
                    continue;
                }
            }
        }

        // Advance sliding window
        if pos + UWINDOW_SIZE < n {
            if seq[pos] == b'T' { t_window = t_window.saturating_sub(1); }
            if pos + UWINDOW_SIZE < n && seq[pos + UWINDOW_SIZE] == b'T' { t_window += 1; }
        }
        pos += 1;
    }

    // Deduplicate overlapping (keep best)
    terminators.sort_by_key(|t| t.start);
    let mut deduped: Vec<Terminator> = Vec::new();
    for t in terminators {
        if let Some(last) = deduped.last() {
            if t.start < last.end {
                if t.confidence > last.confidence {
                    *deduped.last_mut().unwrap() = t;
                }
                continue;
            }
        }
        deduped.push(t);
    }
    deduped
}

/// Find terminators on both strands.
pub fn find_all_terminators(genome: &[u8], rc: &[u8]) -> Vec<Terminator> {
    let glen = genome.len();
    let mut all = find_terminators(genome, true, glen);
    let mut rev = find_terminators(rc, false, glen);
    all.append(&mut rev);
    all.sort_by_key(|t| t.start);
    all
}
