use crate::types::{Gene, MIN_ORF};

#[inline(always)]
pub fn is_start(s: &[u8]) -> bool {
    s == b"ATG" || s == b"GTG" || s == b"TTG"
}

/// Genetic code variants for stop codon recognition.
/// Code 11: Standard bacterial (TAA, TAG, TGA)
/// Code 4:  Mycoplasma/Spiroplasma (TAA, TAG only — TGA = Trp)
#[derive(Clone, Copy, PartialEq)]
pub enum GeneticCode {
    Standard,   // code 11: TAA, TAG, TGA
    Code4,      // code 4:  TAA, TAG (TGA = Trp)
}

#[inline(always)]
pub fn is_stop(s: &[u8]) -> bool {
    s == b"TAA" || s == b"TAG" || s == b"TGA"
}

#[inline(always)]
pub fn is_stop_code(s: &[u8], code: GeneticCode) -> bool {
    match code {
        GeneticCode::Standard => s == b"TAA" || s == b"TAG" || s == b"TGA",
        GeneticCode::Code4 => s == b"TAA" || s == b"TAG",
    }
}

pub fn find_orfs(seq: &[u8], is_plus: bool, glen: usize) -> Vec<Gene> {
    find_orfs_code(seq, is_plus, glen, GeneticCode::Standard)
}

pub fn find_orfs_code(seq: &[u8], is_plus: bool, glen: usize, code: GeneticCode) -> Vec<Gene> {
    let mut orfs = Vec::new();
    let n = seq.len();
    let mut sg = 0u32;

    for frame in 0..3u8 {
        let fr = frame as usize;
        let mut stops = Vec::new();
        let mut i = fr;
        while i + 2 < n {
            if is_stop_code(&seq[i..i+3], code) { stops.push(i); }
            i += 3;
        }

        let mut prev = fr;
        for &sp in &stops {
            let se = sp + 3;
            let mut starts = Vec::new();
            let mut j = prev;
            while j + 2 < n && j < sp {
                if is_start(&seq[j..j+3]) { starts.push(j); }
                j += 3;
            }

            for (k, &ss) in starts.iter().enumerate() {
                let len = se - ss;
                if len < MIN_ORF { continue; }

                let (gs, ge) = if is_plus {
                    (ss + 1, se)
                } else {
                    (glen - se + 1, glen - ss)
                };

                let mut g = Gene::new();
                g.start = gs;
                g.end = ge;
                g.is_plus = is_plus;
                g.seq_start = ss;
                g.seq_end = se;
                g.length = len;
                g.stop_group = sg;
                g.frame = frame;
                g.start_codon.copy_from_slice(&seq[ss..ss+3]);
                g.is_longest = k == 0;
                orfs.push(g);
            }
            sg += 1;
            prev = se;
        }
    }
    orfs
}

/// Auto-detect genetic code by comparing ORF statistics.
///
/// Method: compare longest ORF lengths under standard vs code 4.
/// In code 4 genomes, standard code breaks genes at internal TGA →
/// ORFs are much shorter. If code 4 gives significantly longer ORFs,
/// the genome uses code 4.
pub fn detect_genetic_code(seq: &[u8]) -> GeneticCode {
    let n = seq.len();
    if n < 10000 { return GeneticCode::Standard; }

    // Find top ORF lengths under both codes
    let lens_std = longest_orf_lengths(seq, GeneticCode::Standard, 200);
    let lens_c4 = longest_orf_lengths(seq, GeneticCode::Code4, 200);

    if lens_std.is_empty() || lens_c4.is_empty() {
        return GeneticCode::Standard;
    }

    // Compare median of top 200 longest ORFs
    let median_std = lens_std[lens_std.len() / 2] as f64;
    let median_c4 = lens_c4[lens_c4.len() / 2] as f64;

    // If code 4 gives >40% longer median ORF, the genome uses code 4
    // (in Mycoplasma: TGA breaks ~72% of genes, so standard code ORFs are much shorter)
    let ratio = median_c4 / median_std.max(1.0);
    eprintln!("  Code detection: median ORF standard={}bp, code4={}bp, ratio={:.2}",
        median_std as usize, median_c4 as usize, ratio);

    // Code 4 detection: Mycoplasma/Spiroplasma use TGA as Trp, not stop.
    //
    // Key insight: in code-4 genomes, the LONGEST ORFs found under code-4
    // contain many in-frame TGA codons (used as Trp). In standard-code
    // genomes, code-4 ORFs are just artifacts where TGA happens to be absent.
    //
    // Test: among the longest code-4 ORFs, count in-frame TGA per 1000 codons.
    // Code-4 genome: high TGA density (TGA is a normal Trp codon, ~1-2%)
    // Standard genome high-GC: low TGA density (TGA is rare because AT-rich)
    // Standard genome normal: moderate TGA density (absent by construction of code-4 ORFs)
    // Simple and robust: code 4 organisms (Mycoplasma, Spiroplasma) are
    // ALWAYS low-GC (25-35%). No known high-GC organism uses code 4.
    // The length ratio inflates at high GC because TGA is AT-rich and rare.
    // So: require low GC + significant length ratio.
    let gc = {
        let gc_count = seq.iter().filter(|&&b| b == b'G' || b == b'C').count();
        gc_count as f64 / n as f64
    };
    eprintln!("  Genome GC={:.1}%, ratio={:.2}", gc * 100.0, ratio);

    if ratio > 1.20 && gc < 0.38 {
        GeneticCode::Code4
    } else {
        GeneticCode::Standard
    }
}

/// Count in-frame TGA density inside the longest code-4 ORFs.
/// Returns TGA per codon. In code-4 genomes TGA is Trp (~1-2% of codons).
/// In standard genomes, code-4 ORFs exist because TGA is absent → density ~0.
fn count_inframe_tga_density(seq: &[u8], sorted_lens: &[usize], code: GeneticCode) -> f64 {
    if sorted_lens.is_empty() { return 0.0; }

    // Use top 50 longest ORFs
    let cutoff_idx = sorted_lens.len().saturating_sub(50);
    let min_len = sorted_lens[cutoff_idx];
    let n = seq.len();

    let mut tga_count = 0u64;
    let mut total_codons = 0u64;

    for frame in 0..3usize {
        let mut orf_start: Option<usize> = None;
        let mut i = frame;
        while i + 3 <= n {
            let codon = &seq[i..i + 3];
            if orf_start.is_none() && is_start(codon) {
                orf_start = Some(i);
            }
            if let Some(start) = orf_start {
                if is_stop_code(codon, code) {
                    let len = i - start;
                    if len >= min_len {
                        // Count TGA inside this ORF (excluding the stop itself)
                        let mut j = start;
                        while j < i {
                            if &seq[j..j + 3] == b"TGA" {
                                tga_count += 1;
                            }
                            total_codons += 1;
                            j += 3;
                        }
                    }
                    orf_start = None;
                }
            }
            i += 3;
        }
    }

    if total_codons > 0 { tga_count as f64 / total_codons as f64 } else { 0.0 }
}

/// Count stop codon usage among the longest ORFs (by sorted lengths).
/// Takes the pre-computed sorted lengths to determine a length threshold
/// (top 100 or top 25%), then counts stop types above that threshold.
fn count_stop_usage_top(seq: &[u8], sorted_lens: &[usize]) -> (u32, u32, u32) {
    if sorted_lens.is_empty() {
        return (0, 0, 0);
    }
    // Use the top 100 longest ORFs or top 25%, whichever gives more
    let cutoff_idx = sorted_lens.len().saturating_sub(100).min(sorted_lens.len() * 3 / 4);
    let min_orf_len = sorted_lens[cutoff_idx];
    count_stop_usage(seq, min_orf_len)
}

/// Count which stop codons terminate long ORFs (standard code).
/// Returns (TAA_count, TAG_count, TGA_count).
fn count_stop_usage(seq: &[u8], min_orf_len: usize) -> (u32, u32, u32) {
    let n = seq.len();
    let mut n_taa = 0u32;
    let mut n_tag = 0u32;
    let mut n_tga = 0u32;

    for frame in 0..3usize {
        let mut orf_start: Option<usize> = None;
        let mut i = frame;
        while i + 3 <= n {
            let codon = &seq[i..i+3];
            if orf_start.is_none() && is_start(codon) {
                orf_start = Some(i);
            }
            if let Some(start) = orf_start {
                if is_stop_code(codon, GeneticCode::Standard) {
                    let len = i - start;
                    if len >= min_orf_len {
                        match codon {
                            b"TAA" => n_taa += 1,
                            b"TAG" => n_tag += 1,
                            b"TGA" => n_tga += 1,
                            _ => {}
                        }
                    }
                    orf_start = None;
                }
            }
            i += 3;
        }
    }

    (n_taa, n_tag, n_tga)
}

/// Get sorted lengths of the N longest ORFs under a given genetic code.
fn longest_orf_lengths(seq: &[u8], code: GeneticCode, top_n: usize) -> Vec<usize> {
    let n = seq.len();
    let mut lengths = Vec::new();

    for frame in 0..3usize {
        let mut orf_start: Option<usize> = None;
        let mut i = frame;
        while i + 3 <= n {
            let codon = &seq[i..i+3];
            if orf_start.is_none() && is_start(codon) {
                orf_start = Some(i);
            }
            if let Some(start) = orf_start {
                if is_stop_code(codon, code) {
                    let len = i - start;
                    if len >= 300 {
                        lengths.push(len);
                    }
                    orf_start = None;
                }
            }
            i += 3;
        }
    }

    lengths.sort_unstable_by(|a, b| b.cmp(a)); // descending
    lengths.truncate(top_n);
    lengths.sort_unstable(); // ascending for median
    lengths
}
