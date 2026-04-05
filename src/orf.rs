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

    if ratio > 1.40 {
        GeneticCode::Code4
    } else {
        GeneticCode::Standard
    }
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
