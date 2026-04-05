//! IS element detection and correction (post-annotation module).
//!
//! De novo approach:
//! 1. Cluster predicted proteins by k-mer similarity → repeat families
//! 2. For each family: find longest copy as reference
//! 3. Extend short copies to match reference length
//! 4. Scan gaps for broken/missed copies

use std::collections::HashMap;
use crate::types::Gene;
use crate::io;

const K: usize = 7; // protein k-mer for clustering

// Reuse codon translation
fn translate_codon(a: u8, b: u8, c: u8) -> u8 {
    match (a, b, c) {
        (b'T',b'T',b'T')|(b'T',b'T',b'C') => b'F', (b'T',b'T',b'A')|(b'T',b'T',b'G') => b'L',
        (b'C',b'T',b'T')|(b'C',b'T',b'C')|(b'C',b'T',b'A')|(b'C',b'T',b'G') => b'L',
        (b'A',b'T',b'T')|(b'A',b'T',b'C')|(b'A',b'T',b'A') => b'I', (b'A',b'T',b'G') => b'M',
        (b'G',b'T',b'T')|(b'G',b'T',b'C')|(b'G',b'T',b'A')|(b'G',b'T',b'G') => b'V',
        (b'T',b'C',b'T')|(b'T',b'C',b'C')|(b'T',b'C',b'A')|(b'T',b'C',b'G') => b'S',
        (b'C',b'C',b'T')|(b'C',b'C',b'C')|(b'C',b'C',b'A')|(b'C',b'C',b'G') => b'P',
        (b'A',b'C',b'T')|(b'A',b'C',b'C')|(b'A',b'C',b'A')|(b'A',b'C',b'G') => b'T',
        (b'G',b'C',b'T')|(b'G',b'C',b'C')|(b'G',b'C',b'A')|(b'G',b'C',b'G') => b'A',
        (b'T',b'A',b'T')|(b'T',b'A',b'C') => b'Y',
        (b'C',b'A',b'T')|(b'C',b'A',b'C') => b'H', (b'C',b'A',b'A')|(b'C',b'A',b'G') => b'Q',
        (b'A',b'A',b'T')|(b'A',b'A',b'C') => b'N', (b'A',b'A',b'A')|(b'A',b'A',b'G') => b'K',
        (b'G',b'A',b'T')|(b'G',b'A',b'C') => b'D', (b'G',b'A',b'A')|(b'G',b'A',b'G') => b'E',
        (b'T',b'G',b'T')|(b'T',b'G',b'C') => b'C', (b'T',b'G',b'G') => b'W',
        (b'C',b'G',b'T')|(b'C',b'G',b'C')|(b'C',b'G',b'A')|(b'C',b'G',b'G') => b'R',
        (b'A',b'G',b'T')|(b'A',b'G',b'C') => b'S', (b'A',b'G',b'A')|(b'A',b'G',b'G') => b'R',
        (b'G',b'G',b'T')|(b'G',b'G',b'C')|(b'G',b'G',b'A')|(b'G',b'G',b'G') => b'G',
        _ => b'X',
    }
}

fn translate(dna: &[u8]) -> Vec<u8> {
    let mut prot = Vec::new();
    let mut i = 0;
    while i + 3 <= dna.len() {
        let aa = translate_codon(dna[i], dna[i+1], dna[i+2]);
        if aa == b'X' { break; }
        prot.push(aa);
        i += 3;
    }
    prot
}

fn get_protein(genome: &[u8], gene: &Gene) -> Vec<u8> {
    let dna = if gene.is_plus {
        genome[gene.start.saturating_sub(1)..gene.end.min(genome.len())].to_vec()
    } else {
        io::rev_comp(&genome[gene.start.saturating_sub(1)..gene.end.min(genome.len())])
    };
    translate(&dna)
}

/// Detect Terminal Inverted Repeats (TIR) around a gene region.
/// Extends `flank` bp on each side, searches for inverted repeat 10-35bp, ≥75% identity.
/// Returns (tir_length, is_start, is_end) if found.
fn detect_tir(genome: &[u8], gene_start: usize, gene_end: usize, flank: usize) -> Option<(usize, usize, usize)> {
    let ext_start = gene_start.saturating_sub(flank);
    let ext_end = (gene_end + flank).min(genome.len());
    if ext_end <= ext_start + 30 { return None; }
    let region = &genome[ext_start..ext_end];
    let rlen = region.len();
    let rc_buf = io::rev_comp(region);

    let mut best_score = 0u32;
    let mut best_result: Option<(usize, usize, usize)> = None;

    // Search: left TIR within first (flank+20)bp, right TIR within last (flank+20)bp
    let search_left = (flank + 20).min(rlen / 2);
    let search_right_start = rlen.saturating_sub(flank + 20);

    for tir_len in 10..=30 {
        for lp in 0..search_left.saturating_sub(tir_len) {
            let left = &region[lp..lp + tir_len];
            // Right TIR = reverse complement, search from end
            for rp in search_right_start..rlen.saturating_sub(tir_len) {
                // RC of right tip
                let right_rc_start = rlen - rp - tir_len;
                if right_rc_start + tir_len > rlen { continue; }
                let right_rc = &rc_buf[right_rc_start..right_rc_start + tir_len];

                let matches = left.iter().zip(right_rc.iter()).filter(|(a, b)| a == b).count() as u32;
                if matches >= (tir_len as u32 * 7 / 10) && matches > best_score {
                    best_score = matches;
                    let is_start = ext_start + lp;
                    let is_end = ext_start + rp + tir_len;
                    best_result = Some((tir_len, is_start, is_end));
                }
            }
        }
    }
    best_result
}

/// Build a start context PWM from repeat family copies.
/// Takes the 60bp window around each copy's start, builds position-frequency matrix.
/// Returns log-likelihood PWM if enough copies.
fn build_repeat_start_pwm(genome: &[u8], genes: &[Gene], members: &[usize], window: usize) -> Option<Vec<[f64; 4]>> {
    if members.len() < 3 { return None; }
    let rc_genome = io::rev_comp(genome);
    let mut counts = vec![[0u32; 4]; window * 2];
    let mut n_seqs = 0u32;

    for &mi in members {
        let g = &genes[mi];
        let seq = if g.is_plus { genome } else { &rc_genome };
        let ss = g.seq_start;
        if ss < window || ss + window > seq.len() { continue; }
        let region = &seq[ss - window..ss + window];
        if region.iter().any(|&b| io::nt4(b) > 3) { continue; }
        for (i, &nt) in region.iter().enumerate() {
            counts[i][io::nt4(nt)] += 1;
        }
        n_seqs += 1;
    }
    if n_seqs < 3 { return None; }

    // Background nucleotide frequency
    let mut bg = [0u64; 4];
    let limit = genome.len().min(500000);
    for i in 0..limit {
        let n = io::nt4(genome[i]);
        if n < 4 { bg[n] += 1; }
    }
    let bg_total = bg.iter().sum::<u64>() as f64;

    let mut pwm = Vec::with_capacity(window * 2);
    for i in 0..window * 2 {
        let mut row = [0.0f64; 4];
        for j in 0..4 {
            let freq = (counts[i][j] as f64 + 0.5) / (n_seqs as f64 + 2.0);
            let bg_freq = (bg[j] as f64 + 0.5) / (bg_total + 2.0);
            row[j] = (freq / bg_freq).ln();
        }
        pwm.push(row);
    }
    Some(pwm)
}

/// Score a candidate start position using a repeat-family PWM.
fn score_repeat_start(seq: &[u8], pos: usize, pwm: &[[f64; 4]], window: usize) -> f64 {
    if pos < window || pos + window > seq.len() { return f64::NEG_INFINITY; }
    let region = &seq[pos - window..pos + window];
    let mut score = 0.0;
    for (i, &nt) in region.iter().enumerate() {
        let n = io::nt4(nt);
        if n > 3 { return f64::NEG_INFINITY; }
        if i < pwm.len() { score += pwm[i][n]; }
    }
    score
}

/// Find repeat families, detect TIRs, and correct truncated copies.
/// Returns number of genes corrected.
pub fn correct_repeat_families(genome: &[u8], genes: &mut Vec<Gene>) -> usize {
    if genes.len() < 20 { return 0; }
    let glen = genome.len();

    // Step 1: Get all proteins
    let proteins: Vec<Vec<u8>> = genes.iter().map(|g| get_protein(genome, g)).collect();

    // Step 2: Build k-mer index
    let mut idx: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for (i, prot) in proteins.iter().enumerate() {
        if prot.len() < K { continue; }
        for j in 0..prot.len() - K + 1 {
            idx.entry(&prot[j..j+K]).or_default().push(i);
        }
    }

    // Step 3: Cluster by similarity (>50% shared k-mers)
    let mut cluster_id = vec![usize::MAX; genes.len()];
    let mut n_clusters = 0usize;

    for i in 0..genes.len() {
        if proteins[i].len() < 30 { continue; }
        let total = proteins[i].len() - K + 1;
        if total == 0 { continue; }

        let mut hits = vec![0u32; genes.len()];
        for j in 0..total {
            if let Some(ids) = idx.get(&proteins[i][j..j+K]) {
                for &si in ids {
                    if si != i { hits[si] += 1; }
                }
            }
        }

        // Find similar genes (>50% shared k-mers)
        let similar: Vec<usize> = hits.iter().enumerate()
            .filter(|(si, &c)| *si != i && c as usize > total / 2)
            .map(|(si, _)| si)
            .collect();

        if similar.is_empty() { continue; }

        // Assign cluster
        let my_cid = if cluster_id[i] != usize::MAX {
            cluster_id[i]
        } else {
            let cid = n_clusters;
            n_clusters += 1;
            cluster_id[i] = cid;
            cid
        };

        for &si in &similar {
            if cluster_id[si] == usize::MAX {
                cluster_id[si] = my_cid;
            }
            // Merge clusters if different
            if cluster_id[si] != my_cid {
                let old = cluster_id[si];
                for c in cluster_id.iter_mut() {
                    if *c == old { *c = my_cid; }
                }
            }
        }
    }

    // Step 4: For each cluster ≥3 copies, find longest and extend short ones
    let mut clusters: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, &cid) in cluster_id.iter().enumerate() {
        if cid != usize::MAX {
            clusters.entry(cid).or_default().push(i);
        }
    }

    // Build GLOBAL start-context PWM from ALL repeat family consensus starts.
    // One RNA polymerase → one start context model for all IS elements.
    let pwm_window = 30;
    let all_repeat_members: Vec<usize> = clusters.values()
        .filter(|m| m.len() >= 3)
        .flat_map(|m| m.iter().copied())
        .collect();
    let global_pwm = build_repeat_start_pwm(genome, genes, &all_repeat_members, pwm_window);

    let mut corrected = 0;

    for (_, members) in &clusters {
        if members.len() < 3 { continue; }

        // Safety: only extend copies with the MODAL length (most common).
        // This avoids extending pseudogene fragments or already-correct copies.
        // Find the most common length (±10% tolerance)
        let lengths: Vec<usize> = members.iter().map(|&i| genes[i].length).collect();
        let modal_len = {
            let mut best_count = 0;
            let mut best_len = 0;
            for &l in &lengths {
                let count = lengths.iter().filter(|&&x| {
                    let diff = if x > l { x - l } else { l - x };
                    diff < l / 10 + 10
                }).count();
                if count > best_count { best_count = count; best_len = l; }
            }
            best_len
        };
        // Need ≥3 copies at modal length to be confident
        let modal_count = lengths.iter().filter(|&&l| {
            let diff = if l > modal_len { l - modal_len } else { modal_len - l };
            diff < modal_len / 10 + 10
        }).count();
        if modal_count < 3 { continue; }
        // Only extend copies AT modal length (skip outliers like pseudogene fragments)
        let extend_these: Vec<usize> = members.iter()
            .filter(|&&i| {
                let diff = if genes[i].length > modal_len { genes[i].length - modal_len } else { modal_len - genes[i].length };
                diff < modal_len / 10 + 10
            })
            .copied().collect();

        // TIR-based extension: detect IS boundary, extend gene to fill it
        let rc_genome_buf = io::rev_comp(genome);

        for &mi in &extend_these {
            let g = &genes[mi];
            let g_start_0 = g.start.saturating_sub(1);
            let g_end_0 = g.end;

            let tir = detect_tir(genome, g_start_0, g_end_0, 150);
            if tir.is_none() { continue; }
            let (_tir_len, _is_start, is_end) = tir.unwrap();
            let is_len = is_end - _is_start;
            if is_len <= g.length + 30 { continue; }

            let seq: &[u8] = if g.is_plus { genome } else { &rc_genome_buf };
            let ss = g.seq_start;
            let se = g.seq_end;
            let extension = is_len - g.length;

            // Find farthest upstream start codon within IS boundary
            let mut best_start: Option<usize> = None;
            let mut pos = ss;
            loop {
                if pos < 3 { break; }
                pos -= 3;
                if ss - pos > extension + 90 { break; }
                if pos + 3 > seq.len() { break; }
                let codon = &seq[pos..pos+3];
                if codon == b"TAA" || codon == b"TAG" || codon == b"TGA" { break; }
                if codon == b"ATG" || codon == b"GTG" || codon == b"TTG" {
                    best_start = Some(pos);
                }
            }

            if let Some(new_ss) = best_start {
                let new_len = se - new_ss;
                let actual_ext = new_len - g.length;
                // Only extend if gain ≥90bp (avoid micro-extensions)
                if actual_ext >= 90 && new_len <= is_len + 30 {
                    let g = &mut genes[mi];
                    g.seq_start = new_ss;
                    g.length = new_len;
                    if new_ss + 3 <= seq.len() {
                        g.start_codon = [seq[new_ss], seq[new_ss+1], seq[new_ss+2]];
                    }
                    if g.is_plus {
                        g.start = new_ss + 1;
                    } else {
                        g.end = glen - new_ss;
                    }
                    corrected += 1;
                }
            }
        }
    }

    if corrected > 0 {
        eprintln!("IS/repeat correction: {} genes extended to match family consensus", corrected);
    }

    corrected
}
