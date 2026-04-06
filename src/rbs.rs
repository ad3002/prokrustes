use crate::types::Gene;
use crate::io::nt4;

static SD_MOTIFS: &[(&[u8], f64)] = &[
    (b"TAAGGAGG", 1.00), (b"AAGGAGG", 1.00),
    (b"AGGAGG", 0.98), (b"AGGAG", 0.90), (b"GGAGG", 0.87),
    (b"TAAGGA", 0.80), (b"AAGGA", 0.76), (b"GAGGT", 0.70),
    (b"AGGA", 0.68), (b"GAGG", 0.68), (b"GGAG", 0.60),
    (b"AGG", 0.40), (b"GGA", 0.35), (b"GAG", 0.28), (b"GG", 0.14),
];

pub fn score_rbs_at(seq: &[u8], start_pos: usize) -> f64 {
    if start_pos < 4 { return 0.0; }
    let up_begin = start_pos.saturating_sub(24);
    let up_end = start_pos.saturating_sub(3);
    if up_end <= up_begin { return 0.0; }
    let upstream = &seq[up_begin..up_end];

    let mut best = 0.0;
    for &(pat, base_score) in SD_MOTIFS {
        if pat.len() > upstream.len() { continue; }
        for i in 0..=upstream.len() - pat.len() {
            if &upstream[i..i+pat.len()] == pat {
                let motif_end = up_begin + i + pat.len();
                let spacing = start_pos as isize - motif_end as isize;
                let sf = if spacing >= 5 && spacing <= 10 { 1.0 }
                    else if spacing >= 4 && spacing <= 12 { 0.85 }
                    else if spacing >= 3 && spacing <= 14 { 0.68 }
                    else if spacing >= 1 && spacing <= 18 { 0.42 }
                    else { 0.12 };
                let s = base_score * sf;
                if s > best { best = s; }
            }
        }
    }
    best
}

/// Build PWM from a specific set of genes (for iterative training).
pub fn build_rbs_pwm_from_genes(seq: &[u8], genes: &[&Gene], width: usize) -> Option<Vec<[f64; 4]>> {
    if genes.len() < 60 { return None; }

    let mut counts = vec![[0u32; 4]; width];
    let mut bg = [0u64; 4];
    let mut bg_total = 0u64;
    let mut n_seqs = 0u32;

    for o in genes {
        if o.seq_start < width + 3 { continue; }
        let region = &seq[o.seq_start - width..o.seq_start];
        if region.iter().any(|&b| nt4(b) > 3) { continue; }
        for (i, &nt) in region.iter().enumerate() {
            counts[i][nt4(nt)] += 1;
        }
        n_seqs += 1;
    }
    if n_seqs < 40 { return None; }

    let limit = seq.len().min(500000);
    for i in 0..limit {
        let n = nt4(seq[i]);
        if n < 4 { bg[n] += 1; bg_total += 1; }
    }
    if bg_total == 0 { return None; }

    let mut pwm = Vec::with_capacity(width);
    for i in 0..width {
        let mut row = [0.0f64; 4];
        for j in 0..4 {
            let freq = (counts[i][j] as f64 + 0.5) / (n_seqs as f64 + 2.0);
            let bg_freq = (bg[j] as f64 + 0.5) / (bg_total as f64 + 2.0);
            row[j] = (freq / bg_freq).ln();
        }
        pwm.push(row);
    }
    Some(pwm)
}

/// Build initial PWM from confident ORFs (first pass).
pub fn build_rbs_pwm(seq: &[u8], orfs: &[Gene]) -> Option<Vec<[f64; 4]>> {
    let conf: Vec<&Gene> = orfs.iter().filter(|o| {
        o.score > 0.46 && o.length >= 350 && o.is_atg() && o.is_longest
    }).collect();
    if conf.len() < 80 { return None; }

    let width = 25usize; // expanded from 20 to 25 (Prodigal uses 25)
    let mut counts = vec![[0u32; 4]; width];
    let mut bg = [0u64; 4];
    let mut bg_total = 0u64;
    let mut n_seqs = 0u32;

    for o in &conf {
        if o.seq_start < width + 3 { continue; }
        let region = &seq[o.seq_start - width..o.seq_start];
        if region.iter().any(|&b| nt4(b) > 3) { continue; }
        for (i, &nt) in region.iter().enumerate() {
            counts[i][nt4(nt)] += 1;
        }
        n_seqs += 1;
    }
    if n_seqs < 60 { return None; }

    let limit = seq.len().min(500000);
    for i in 0..limit {
        let n = nt4(seq[i]);
        if n < 4 { bg[n] += 1; bg_total += 1; }
    }
    if bg_total == 0 { return None; }

    let mut pwm = Vec::with_capacity(width);
    for i in 0..width {
        let mut row = [0.0f64; 4];
        for j in 0..4 {
            let freq = (counts[i][j] as f64 + 0.5) / (n_seqs as f64 + 2.0);
            let bg_freq = (bg[j] as f64 + 0.5) / (bg_total as f64 + 2.0);
            row[j] = (freq / bg_freq).ln();
        }
        pwm.push(row);
    }
    Some(pwm)
}

pub fn score_rbs_pwm(seq: &[u8], start_pos: usize, pwm: &Option<Vec<[f64; 4]>>) -> f64 {
    let pwm = match pwm { Some(p) => p, None => return 0.0 };
    let w = pwm.len();
    if start_pos < w + 3 { return 0.0; }
    let region = &seq[start_pos - w..start_pos];
    if region.iter().any(|&b| nt4(b) > 3) { return 0.0; }
    let mut total = 0.0;
    for (i, &nt) in region.iter().enumerate() {
        total += pwm[i][nt4(nt)];
    }
    total
}

/// Detect leaderless genes via canonical promoter motifs.
///
/// Bacterial promoter = -35 box (TTGACA, ~35bp upstream) + -10 box (TATAAT, ~10bp upstream).
/// Leaderless genes have promoter immediately upstream — transcription starts at/near ATG.
/// SD-dependent genes have SD 5-13bp upstream instead.
///
/// Scoring:
///   -10 box alone (5/6 match): 0.5
///   -35 box alone (5/6 match): 0.3
///   Both -10 and -35 in correct spacing (15-21bp apart): 1.0
///   Partial matches scaled proportionally.
pub fn score_leaderless(seq: &[u8], start_pos: usize) -> f64 {
    if start_pos < 15 { return 0.0; }

    // Search for -10 box (TATAAT) at 3-14bp upstream of start
    let minus10_score = {
        let rs = start_pos.saturating_sub(14);
        let re = start_pos.saturating_sub(3);
        if re <= rs + 5 { 0.0 }
        else {
            let region = &seq[rs..re];
            let consensus = b"TATAAT";
            let mut best = 0.0f64;
            let mut best_pos = 0usize;
            for i in 0..region.len().saturating_sub(5) {
                let m = region[i..i+6].iter().zip(consensus).filter(|(a,b)| a == b).count();
                let s = m as f64 / 6.0;
                if s > best { best = s; best_pos = rs + i; }
            }
            if best >= 4.0 / 6.0 { best } else { 0.0 }
        }
    };

    // Search for -35 box (TTGACA) at 25-45bp upstream of start
    let minus35_score = if start_pos >= 25 {
        let rs = start_pos.saturating_sub(45);
        let re = start_pos.saturating_sub(25);
        if re <= rs + 5 { 0.0 }
        else {
            let region = &seq[rs..re.min(seq.len())];
            let consensus = b"TTGACA";
            let mut best = 0.0f64;
            for i in 0..region.len().saturating_sub(5) {
                let m = region[i..i+6].iter().zip(consensus).filter(|(a,b)| a == b).count();
                let s = m as f64 / 6.0;
                if s > best { best = s; }
            }
            if best >= 4.0 / 6.0 { best } else { 0.0 }
        }
    } else {
        0.0
    };

    // Combined score
    if minus10_score > 0.0 && minus35_score > 0.0 {
        // Both boxes found — strong promoter signal
        (minus10_score * 0.6 + minus35_score * 0.4).min(1.0)
    } else if minus10_score >= 5.0 / 6.0 {
        // Only -10 box (original behavior, but stricter threshold)
        minus10_score * 0.5
    } else if minus35_score >= 5.0 / 6.0 {
        minus35_score * 0.3
    } else {
        0.0
    }
}

/// Score nucleotide context around start codon (cycle 25 idea).
/// Positions -3,-1 prefer A/T, +4 prefers G in E. coli.
pub fn score_start_context(seq: &[u8], start_pos: usize) -> f64 {
    let mut s = 0.0;
    if start_pos >= 3 {
        if seq[start_pos-3] == b'A' || seq[start_pos-3] == b'T' { s += 0.15; }
        if seq[start_pos-1] == b'A' || seq[start_pos-1] == b'T' { s += 0.08; }
        if seq[start_pos-2] == b'A' { s += 0.05; }
    }
    if start_pos + 4 < seq.len() {
        if seq[start_pos+3] == b'G' { s += 0.10; }
        else if seq[start_pos+3] == b'A' { s += 0.05; }
    }
    s
}

pub fn compute_upstream_at(seq: &[u8], start_pos: usize) -> f64 {
    let begin = start_pos.saturating_sub(40);
    if begin >= start_pos || start_pos - begin < 5 { return 0.5; }
    let region = &seq[begin..start_pos];
    let at = region.iter().filter(|&&c| c == b'A' || c == b'T').count();
    at as f64 / region.len() as f64
}

/// Trained upstream composition model.
/// Stores per-position nucleotide preferences at 32 positions relative to start.
/// Positions: -1, -2 (start context) and -15 to -44 (upstream, skipping RBS zone).
pub struct UpstreamModel {
    /// weights[pos_index][nt] where nt: 0=A, 1=C, 2=G, 3=T
    pub weights: [[f64; 4]; 32],
}

impl UpstreamModel {
    /// Train upstream model from confident genes.
    pub fn train(seq: &[u8], genes: &[Gene]) -> Option<Self> {
        if genes.len() < 100 { return None; }

        // Count nucleotide frequencies at each position
        let mut counts = [[0u32; 4]; 32];
        let mut totals = [0u32; 32];

        // Background counts (genome-wide)
        let mut bg = [0u32; 4];
        let mut bg_total = 0u32;
        for &b in seq {
            let n = nt4(b);
            if n < 4 { bg[n] += 1; bg_total += 1; }
        }

        for gene in genes {
            if gene.score < 0.45 || gene.length < 300 { continue; }
            let sp = gene.seq_start;

            // Positions -1, -2 (indices 0, 1)
            for offset in 1..=2usize {
                if sp >= offset {
                    let n = nt4(seq[sp - offset]);
                    if n < 4 {
                        counts[offset - 1][n] += 1;
                        totals[offset - 1] += 1;
                    }
                }
            }

            // Positions -15 to -44 (indices 2..32)
            for i in 0..30usize {
                let pos = 15 + i; // distance from start
                if sp >= pos {
                    let n = nt4(seq[sp - pos]);
                    if n < 4 {
                        counts[i + 2][n] += 1;
                        totals[i + 2] += 1;
                    }
                }
            }
        }

        // Compute log-likelihood ratios vs background
        let mut weights = [[0.0f64; 4]; 32];
        let bg_f = bg_total as f64;

        for p in 0..32 {
            if totals[p] < 50 { continue; }
            let t_f = totals[p] as f64;
            for nt in 0..4 {
                let p_gene = (counts[p][nt] as f64 + 1.0) / (t_f + 4.0);
                let p_bg = (bg[nt] as f64 + 1.0) / (bg_f + 4.0);
                weights[p][nt] = (p_gene / p_bg).ln();
            }
        }

        Some(UpstreamModel { weights })
    }

    /// Score a single start position using the upstream model.
    pub fn score(&self, seq: &[u8], start_pos: usize) -> f64 {
        let mut total = 0.0f64;
        let mut n = 0u32;

        // Positions -1, -2
        for offset in 1..=2usize {
            if start_pos >= offset {
                let nt = nt4(seq[start_pos - offset]);
                if nt < 4 {
                    total += self.weights[offset - 1][nt];
                    n += 1;
                }
            }
        }

        // Positions -15 to -44
        for i in 0..30usize {
            let pos = 15 + i;
            if start_pos >= pos {
                let nt = nt4(seq[start_pos - pos]);
                if nt < 4 {
                    total += self.weights[i + 2][nt];
                    n += 1;
                }
            }
        }

        if n > 0 { total / n as f64 } else { 0.0 }
    }
}
