use crate::types::{Gene, HexModel, MonoModel, N_HEX, N_TRI};
use crate::io::{hex_enc, tri_enc, rev_comp};

/// Build initial hexamer table from long ORFs (does NOT require is_longest).
pub fn train_hex_initial(seq: &[u8], orfs: &[Gene], min_len: usize, require_atg: bool) -> Option<HexModel> {
    let training: Vec<&Gene> = orfs.iter().filter(|o| {
        o.length >= min_len && (!require_atg || o.is_atg())
    }).collect();

    if training.is_empty() { return None; }

    let mut cod = [0u32; N_HEX];
    let mut ct = 0u64;
    let mut ncod = [0u32; N_HEX];
    let mut nt = 0u64;

    for orf in &training {
        let s = &seq[orf.seq_start..orf.seq_end];
        // In-frame hexamers -> coding
        let mut i = 0;
        while i + 5 < s.len() {
            if let Some(idx) = hex_enc(&s[i..i+6]) {
                cod[idx] += 1;
                ct += 1;
            }
            i += 3;
        }
        // Out-of-frame -> non-coding
        for shift in [1usize, 2] {
            let mut i = shift;
            while i + 5 < s.len() {
                if let Some(idx) = hex_enc(&s[i..i+6]) {
                    ncod[idx] += 1;
                    nt += 1;
                }
                i += 3;
            }
        }
        // Reverse complement -> non-coding
        let rcs = rev_comp(s);
        for fo in 0..3 {
            let mut i = fo;
            while i + 5 < rcs.len() {
                if let Some(idx) = hex_enc(&rcs[i..i+6]) {
                    ncod[idx] += 1;
                    nt += 1;
                }
                i += 3;
            }
        }
    }

    if ct < 500 || nt < 500 { return None; }
    Some(build_hex_table(&cod, ct, &ncod, nt))
}

/// Build hexamer table from pre-filtered ORFs (caller does filtering).
pub fn train_hex_from_set(seq: &[u8], orfs: &[Gene], min_len: usize) -> Option<HexModel> {
    let training: Vec<&Gene> = orfs.iter().filter(|o| o.length >= min_len).collect();

    if training.is_empty() { return None; }

    let mut cod = [0u32; N_HEX];
    let mut ct = 0u64;
    let mut ncod = [0u32; N_HEX];
    let mut nt = 0u64;

    for orf in &training {
        let s = &seq[orf.seq_start..orf.seq_end];
        let mut i = 0;
        while i + 5 < s.len() {
            if let Some(idx) = hex_enc(&s[i..i+6]) {
                cod[idx] += 1; ct += 1;
            }
            i += 3;
        }
        for shift in [1usize, 2] {
            let mut i = shift;
            while i + 5 < s.len() {
                if let Some(idx) = hex_enc(&s[i..i+6]) {
                    ncod[idx] += 1; nt += 1;
                }
                i += 3;
            }
        }
        let rcs = rev_comp(s);
        for fo in 0..3 {
            let mut i = fo;
            while i + 5 < rcs.len() {
                if let Some(idx) = hex_enc(&rcs[i..i+6]) {
                    ncod[idx] += 1; nt += 1;
                }
                i += 3;
            }
        }
    }

    if ct < 500 || nt < 500 { return None; }
    Some(build_hex_table(&cod, ct, &ncod, nt))
}

pub fn build_hex_table(cod: &[u32; N_HEX], ct: u64, ncod: &[u32; N_HEX], nt: u64) -> HexModel {
    let mut model = [0.0f64; N_HEX];
    for i in 0..N_HEX {
        let cf = (cod[i] as f64 + 1.0) / (ct as f64 + N_HEX as f64);
        let nf = (ncod[i] as f64 + 1.0) / (nt as f64 + N_HEX as f64);
        model[i] = (cf / nf).ln();
    }
    model
}

/// Cached hex indices for one ORF: (hex_index, frame) pairs.
pub struct HexCacheEntry {
    pub indices: Vec<(u16, u8)>, // (hex_index, frame)
}

/// Pre-computed hex indices for all ORFs on one strand.
pub struct HexCache {
    pub entries: Vec<HexCacheEntry>,
}

impl HexCache {
    /// Build cache: scan each ORF once, store hex indices + frames.
    pub fn build(seq: &[u8], orfs: &[Gene]) -> Self {
        let entries = orfs.iter().map(|orf| {
            let s = &seq[orf.seq_start..orf.seq_end];
            let ns = s.len();
            if ns < 6 {
                return HexCacheEntry { indices: Vec::new() };
            }
            let mut indices = Vec::with_capacity(ns - 5);
            for i in 0..ns.saturating_sub(5) {
                if let Some(idx) = hex_enc(&s[i..i+6]) {
                    indices.push((idx as u16, (i % 3) as u8));
                }
            }
            HexCacheEntry { indices }
        }).collect();
        HexCache { entries }
    }

    /// Rescore all ORFs using cached indices + new model. No sequence access needed.
    pub fn rescore(&self, orfs: &mut [Gene], model: &HexModel) {
        for (orf, entry) in orfs.iter_mut().zip(self.entries.iter()) {
            if entry.indices.is_empty() {
                orf.hex_avg = 0.0; orf.hex_total = 0.0;
                orf.frame_bias = 0.0; orf.hex_cov = 0.5;
                continue;
            }

            let mut scores = [0.0f64; 3];
            let mut counts = [0u32; 3];
            let mut pos_count = 0u32;

            for &(idx, frame) in &entry.indices {
                let f = frame as usize;
                scores[f] += model[idx as usize];
                counts[f] += 1;
                if f == 0 && model[idx as usize] > 0.0 { pos_count += 1; }
            }

            let n0 = counts[0];
            orf.hex_total = scores[0];
            orf.hex_avg = if n0 > 0 { scores[0] / n0 as f64 } else { 0.0 };
            orf.hex_cov = if n0 > 0 { pos_count as f64 / n0 as f64 } else { 0.5 };

            let avgs: [f64; 3] = [
                if counts[0] > 0 { scores[0] / counts[0] as f64 } else { 0.0 },
                if counts[1] > 0 { scores[1] / counts[1] as f64 } else { 0.0 },
                if counts[2] > 0 { scores[2] / counts[2] as f64 } else { 0.0 },
            ];
            orf.frame_bias = avgs[0] - avgs[1].max(avgs[2]);
        }
    }
}

pub fn score_hex_all(seq: &[u8], orfs: &mut [Gene], model: &HexModel) {
    for orf in orfs.iter_mut() {
        let s = &seq[orf.seq_start..orf.seq_end];
        let ns = s.len();
        if ns < 6 {
            orf.hex_avg = 0.0; orf.hex_total = 0.0;
            orf.frame_bias = 0.0; orf.hex_cov = 0.5;
            continue;
        }

        let mut scores = [0.0f64; 3];
        let mut counts = [0u32; 3];
        let mut pos_count = 0u32;

        for i in 0..ns.saturating_sub(5) {
            if let Some(idx) = hex_enc(&s[i..i+6]) {
                let f = i % 3;
                scores[f] += model[idx];
                counts[f] += 1;
                if f == 0 && model[idx] > 0.0 { pos_count += 1; }
            }
        }

        let n0 = counts[0];
        orf.hex_total = scores[0];
        orf.hex_avg = if n0 > 0 { scores[0] / n0 as f64 } else { 0.0 };
        orf.hex_cov = if n0 > 0 { pos_count as f64 / n0 as f64 } else { 0.5 };

        let avgs: [f64; 3] = [
            if counts[0] > 0 { scores[0] / counts[0] as f64 } else { 0.0 },
            if counts[1] > 0 { scores[1] / counts[1] as f64 } else { 0.0 },
            if counts[2] > 0 { scores[2] / counts[2] as f64 } else { 0.0 },
        ];
        orf.frame_bias = avgs[0] - avgs[1].max(avgs[2]);
    }
}

pub fn blend_hex(a: &HexModel, b: &HexModel, w: f64) -> HexModel {
    let mut r = [0.0f64; N_HEX];
    for i in 0..N_HEX {
        r[i] = (1.0 - w) * a[i] + w * b[i];
    }
    r
}

pub fn merge_hex(a: &Option<HexModel>, b: &Option<HexModel>) -> Option<HexModel> {
    match (a, b) {
        (Some(a), Some(b)) => Some(blend_hex(a, b, 0.5)),
        (Some(a), None) => Some(*a),
        (None, Some(b)) => Some(*b),
        (None, None) => None,
    }
}

pub fn train_intergenic_hex(seq: &[u8], coding_orfs: &[Gene], intergenic: &[Vec<u8>]) -> Option<HexModel> {
    if coding_orfs.is_empty() || intergenic.is_empty() { return None; }

    let mut cod = [0u32; N_HEX];
    let mut ct = 0u64;
    for orf in coding_orfs {
        let s = &seq[orf.seq_start..orf.seq_end];
        let mut i = 0;
        while i + 5 < s.len() {
            if let Some(idx) = hex_enc(&s[i..i+6]) {
                cod[idx] += 1; ct += 1;
            }
            i += 3;
        }
    }
    if ct == 0 { return None; }

    let mut ncod = [0u32; N_HEX];
    let mut nt = 0u64;
    for s in intergenic {
        for fo in 0..3 {
            let mut i = fo;
            while i + 5 < s.len() {
                if let Some(idx) = hex_enc(&s[i..i+6]) {
                    ncod[idx] += 1; nt += 1;
                }
                i += 3;
            }
        }
    }
    if nt < 500 { return None; }
    Some(build_hex_table(&cod, ct, &ncod, nt))
}

// --- Monocodon model ---

pub fn train_mono(seq: &[u8], orfs: &[Gene], min_len: usize) -> Option<MonoModel> {
    let training: Vec<&Gene> = orfs.iter().filter(|o| {
        o.length >= min_len && o.is_atg() && o.is_longest
    }).collect();

    if training.len() < 50 { return None; }

    let mut cod = [0u32; N_TRI];
    let mut ct = 0u64;
    let mut ncod = [0u32; N_TRI];
    let mut nt = 0u64;

    for orf in &training {
        let s = &seq[orf.seq_start..orf.seq_end];
        let mut i = 0;
        while i + 2 < s.len() {
            if let Some(idx) = tri_enc(&s[i..i+3]) {
                cod[idx] += 1;
                ct += 1;
            }
            i += 3;
        }
    }

    for orf in &training {
        let s = &seq[orf.seq_start..orf.seq_end];
        for shift in [1usize, 2] {
            let mut i = shift;
            while i + 2 < s.len() {
                if let Some(idx) = tri_enc(&s[i..i+3]) {
                    ncod[idx] += 1;
                    nt += 1;
                }
                i += 3;
            }
        }
    }

    if ct == 0 || nt == 0 { return None; }

    let mut table = [0.0f64; N_TRI];
    for i in 0..N_TRI {
        let cf = (cod[i] as f64 + 1.0) / (ct as f64 + N_TRI as f64);
        let nf = (ncod[i] as f64 + 1.0) / (nt as f64 + N_TRI as f64);
        table[i] = (cf / nf).ln();
    }
    Some(table)
}

pub fn train_mono_from_set(seq: &[u8], orfs: &[Gene], min_len: usize) -> Option<MonoModel> {
    let training: Vec<&Gene> = orfs.iter().filter(|o| o.length >= min_len).collect();
    if training.len() < 50 { return None; }

    let mut cod = [0u32; N_TRI];
    let mut ct = 0u64;
    let mut ncod = [0u32; N_TRI];
    let mut nt = 0u64;

    for orf in &training {
        let s = &seq[orf.seq_start..orf.seq_end];
        let mut i = 0;
        while i + 2 < s.len() {
            if let Some(idx) = tri_enc(&s[i..i+3]) {
                cod[idx] += 1; ct += 1;
            }
            i += 3;
        }
    }
    for orf in &training {
        let s = &seq[orf.seq_start..orf.seq_end];
        for shift in [1usize, 2] {
            let mut i = shift;
            while i + 2 < s.len() {
                if let Some(idx) = tri_enc(&s[i..i+3]) {
                    ncod[idx] += 1; nt += 1;
                }
                i += 3;
            }
        }
    }
    if ct == 0 || nt == 0 { return None; }

    let mut table = [0.0f64; N_TRI];
    for i in 0..N_TRI {
        let cf = (cod[i] as f64 + 1.0) / (ct as f64 + N_TRI as f64);
        let nf = (ncod[i] as f64 + 1.0) / (nt as f64 + N_TRI as f64);
        table[i] = (cf / nf).ln();
    }
    Some(table)
}

pub fn score_mono_all(seq: &[u8], orfs: &mut [Gene], model: &MonoModel) {
    for orf in orfs.iter_mut() {
        let s = &seq[orf.seq_start..orf.seq_end];
        let mut total = 0.0;
        let mut n = 0u32;
        let mut i = 0;
        while i + 2 < s.len() {
            if let Some(idx) = tri_enc(&s[i..i+3]) {
                total += model[idx];
                n += 1;
            }
            i += 3;
        }
        orf.mono = if n > 0 { total / n as f64 } else { 0.0 };
    }
}

pub fn merge_mono(a: &Option<MonoModel>, b: &Option<MonoModel>) -> Option<MonoModel> {
    match (a, b) {
        (Some(a), Some(b)) => {
            let mut r = [0.0f64; N_TRI];
            for i in 0..N_TRI { r[i] = (a[i] + b[i]) / 2.0; }
            Some(r)
        }
        (Some(a), None) => Some(*a),
        (None, Some(b)) => Some(*b),
        (None, None) => None,
    }
}

pub fn blend_mono(a: &MonoModel, b: &MonoModel, w: f64) -> MonoModel {
    let mut r = [0.0f64; N_TRI];
    for i in 0..N_TRI {
        r[i] = (1.0 - w) * a[i] + w * b[i];
    }
    r
}

/// Upstream trimer model (cycle 24 idea): 3-mer log-likelihood for region
/// 25bp upstream of start codon vs distant upstream (100-200bp).
/// Captures start-site specific trinucleotide patterns (e.g. SD-like motifs).
pub fn train_upstream_trimer(seq: &[u8], orfs: &[Gene], min_score: f64) -> Option<MonoModel> {
    let conf: Vec<&Gene> = orfs.iter()
        .filter(|o| o.score > min_score && o.length >= 400 && o.is_atg() && o.is_longest)
        .collect();
    if conf.len() < 100 { return None; }

    let mut pos_counts = [0u32; N_TRI];
    let mut pos_total = 0u32;
    let mut neg_counts = [0u32; N_TRI];
    let mut neg_total = 0u32;

    for o in &conf {
        let sp = o.seq_start;
        if sp < 30 { continue; }
        // Positive: 25bp just upstream of start
        for i in (sp.saturating_sub(25))..sp.saturating_sub(2) {
            if i + 3 <= seq.len() {
                if let Some(idx) = super::io::tri_enc(&seq[i..i+3]) {
                    pos_counts[idx] += 1;
                    pos_total += 1;
                }
            }
        }
        // Negative: distant upstream (100-200bp before start)
        let neg_start = sp.saturating_sub(200);
        let neg_end = sp.saturating_sub(100);
        for i in neg_start..neg_end {
            if i + 3 <= seq.len() {
                if let Some(idx) = super::io::tri_enc(&seq[i..i+3]) {
                    neg_counts[idx] += 1;
                    neg_total += 1;
                }
            }
        }
    }

    if pos_total < 500 || neg_total < 500 { return None; }

    let mut model = [0.0f64; N_TRI];
    for i in 0..N_TRI {
        let pf = (pos_counts[i] as f64 + 0.5) / (pos_total as f64 + 32.0);
        let nf = (neg_counts[i] as f64 + 0.5) / (neg_total as f64 + 32.0);
        model[i] = (pf / nf).ln();
    }
    Some(model)
}

/// Score upstream region with trimer model.
pub fn score_upstream_trimer(seq: &[u8], start_pos: usize, model: &MonoModel) -> f64 {
    if start_pos < 25 { return 0.0; }
    let mut total = 0.0;
    let mut n = 0u32;
    for i in (start_pos.saturating_sub(25))..start_pos.saturating_sub(2) {
        if i + 3 <= seq.len() {
            if let Some(idx) = super::io::tri_enc(&seq[i..i+3]) {
                total += model[idx];
                n += 1;
            }
        }
    }
    if n > 0 { total / n as f64 } else { 0.0 }
}

/// Train atypical (AT-rich) hexamer model (GeneMarkS-2 idea).
/// ~15% of E. coli genes come from horizontal gene transfer and have
/// GC content much lower than genome average (~0.508). Standard hexamer
/// model misses them because it's trained on typical GC-rich genes.
///
/// Strategy: find ORFs with GC < genome_gc - 0.06 that are still long enough
/// to be real genes. Train a separate hexamer model on them.
pub fn train_atypical_hex(seq: &[u8], orfs: &[Gene], genome_gc: f64) -> Option<HexModel> {
    let gc_cutoff = genome_gc - 0.06;

    // Find AT-rich ORFs that are likely real genes (long, ATG start)
    let atypical: Vec<&Gene> = orfs.iter().filter(|o| {
        o.gc3 < gc_cutoff && o.length >= 300 && o.is_atg() && o.is_longest
    }).collect();

    if atypical.len() < 30 { return None; }

    // Count hexamers in atypical coding regions
    let mut cod = [0u32; N_HEX];
    let mut cod_total = 0u64;
    for orf in &atypical {
        let s = &seq[orf.seq_start..orf.seq_end];
        for i in (0..s.len().saturating_sub(5)).step_by(3) {
            if let Some(idx) = hex_enc(&s[i..i+6]) {
                cod[idx] += 1;
                cod_total += 1;
            }
        }
    }
    if cod_total < 5000 { return None; }

    // Background: non-coding regions (all frames shuffled)
    let mut bg = [0u32; N_HEX];
    let mut bg_total = 0u64;
    for offset in 0..3 {
        for i in (offset..seq.len().saturating_sub(5)).step_by(3) {
            if let Some(idx) = hex_enc(&seq[i..i+6]) {
                bg[idx] += 1;
                bg_total += 1;
            }
        }
    }
    if bg_total == 0 { return None; }

    Some(build_hex_table(&cod, cod_total, &bg, bg_total))
}

/// Score ORFs with atypical model. Only affects low-GC ORFs.
/// Returns improvement over standard model for AT-rich ORFs.
pub fn score_atypical(seq: &[u8], orfs: &mut [Gene], model: &HexModel, genome_gc: f64) {
    let gc_cutoff = genome_gc - 0.04;
    for orf in orfs.iter_mut() {
        if orf.gc3 >= gc_cutoff { continue; } // only score AT-rich ORFs
        let s = &seq[orf.seq_start..orf.seq_end];
        if s.len() < 6 { continue; }
        let mut total = 0.0f64;
        let mut n = 0u32;
        for i in (0..s.len().saturating_sub(5)).step_by(3) {
            if let Some(idx) = hex_enc(&s[i..i+6]) {
                total += model[idx];
                n += 1;
            }
        }
        let atypical_score = if n > 0 { total / n as f64 } else { 0.0 };
        // If atypical model scores this ORF better than standard, use atypical score
        if atypical_score > orf.hex_avg {
            orf.hex_avg = atypical_score;
            orf.score = (orf.score + 0.06).min(1.0);
        }
    }
}

// --- Dicodon (codon-pair) scoring ---
// A dicodon is 2 adjacent codons read in-frame (6bp at codon boundaries).
// Same size as hexamer (4096 entries) but captures codon-pair bias rather
// than sliding-window hexamer frequencies.

/// Train dicodon model from confident genes.
pub fn train_dicodon(seq: &[u8], orfs: &[Gene], min_len: usize) -> Option<HexModel> {
    let training: Vec<&Gene> = orfs.iter()
        .filter(|o| o.length >= min_len && o.is_longest)
        .collect();
    if training.len() < 50 { return None; }

    let mut cod = [0u32; N_HEX]; // coding dicodon counts
    let mut ct = 0u64;
    let mut ncod = [0u32; N_HEX]; // noncoding dicodon counts (from intergenic)
    let mut nt = 0u64;

    // Coding: walk in-frame through training genes
    for orf in &training {
        let s = &seq[orf.seq_start..orf.seq_end];
        let mut i = 0;
        while i + 6 <= s.len() {
            if let Some(idx) = hex_enc(&s[i..i+6]) {
                cod[idx] += 1;
                ct += 1;
            }
            i += 6; // step by 2 codons (6bp) — always in-frame
        }
    }

    // Noncoding: sample from intergenic regions
    // Use regions between training genes
    let mut gene_intervals: Vec<(usize, usize)> = training.iter()
        .map(|o| (o.seq_start, o.seq_end)).collect();
    gene_intervals.sort();

    let mut prev_end = 0;
    for &(start, end) in &gene_intervals {
        if start > prev_end + 50 {
            let ig = &seq[prev_end..start];
            // Sample dicodons from all 3 frames in intergenic
            for frame in 0..3 {
                let mut i = frame;
                while i + 6 <= ig.len() {
                    if let Some(idx) = hex_enc(&ig[i..i+6]) {
                        ncod[idx] += 1;
                        nt += 1;
                    }
                    i += 6;
                }
            }
        }
        prev_end = prev_end.max(end);
    }

    if ct < 1000 || nt < 1000 { return None; }

    // Log-likelihood ratio
    let mut model = [0.0f64; N_HEX];
    let ct_f = ct as f64;
    let nt_f = nt as f64;
    for i in 0..N_HEX {
        let p_cod = (cod[i] as f64 + 1.0) / (ct_f + N_HEX as f64);
        let p_ncod = (ncod[i] as f64 + 1.0) / (nt_f + N_HEX as f64);
        model[i] = (p_cod / p_ncod).ln();
    }
    Some(model)
}

/// Score all ORFs with dicodon model (in-frame only).
pub fn score_dicodon_all(seq: &[u8], orfs: &mut [Gene], model: &HexModel) {
    for orf in orfs.iter_mut() {
        let s = &seq[orf.seq_start..orf.seq_end];
        let mut total = 0.0f64;
        let mut n = 0u32;
        let mut i = 0;
        while i + 6 <= s.len() {
            if let Some(idx) = hex_enc(&s[i..i+6]) {
                total += model[idx];
                n += 1;
            }
            i += 6;
        }
        orf.dicodon = if n > 0 { total / n as f64 } else { 0.0 };
    }
}
