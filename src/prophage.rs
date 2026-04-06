//! Prophage detection and routing module.
//!
//! Detects prophage islands de novo using:
//! 1. k-mer composition divergence (KL divergence vs genome background)
//! 2. Coding density anomalies (phages pack tighter)
//! 3. GC content shifts
//!
//! Returns per-window classification: Host / Prophage / Cryptic
//! Used by the main pipeline to route regions to different scoring priors.

use crate::io::nt4;

const WINDOW: usize = 8000;   // 8kb sliding window
const STEP: usize = 2000;     // 2kb step
const KMER_K: usize = 4;      // tetranucleotide frequencies
const N_KMERS: usize = 256;   // 4^4

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum RegionType {
    Host,
    Prophage,
    Cryptic, // degraded prophage
}

/// Per-window feature vector for region classification.
pub struct WindowFeatures {
    pub start: usize,
    pub end: usize,
    pub kl_divergence: f64,   // KL(window || genome) — higher = more alien
    pub gc_shift: f64,        // |local_gc - genome_gc|
    pub coding_density: f64,  // fraction of window covered by ORFs
    pub mean_orf_length: f64, // average ORF length in window
    pub region_type: RegionType,
}

/// Compute k-mer frequency profile for a sequence.
fn kmer_profile(seq: &[u8]) -> [f64; N_KMERS] {
    let mut counts = [0u32; N_KMERS];
    let mut total = 0u32;

    for i in 0..seq.len().saturating_sub(KMER_K - 1) {
        let mut idx = 0usize;
        let mut valid = true;
        for j in 0..KMER_K {
            let n = nt4(seq[i + j]);
            if n > 3 { valid = false; break; }
            idx = idx * 4 + n;
        }
        if valid {
            counts[idx] += 1;
            total += 1;
        }
    }

    let mut profile = [0.0f64; N_KMERS];
    let t = total as f64;
    if t > 0.0 {
        for i in 0..N_KMERS {
            profile[i] = (counts[i] as f64 + 0.5) / (t + N_KMERS as f64 * 0.5);
        }
    } else {
        for i in 0..N_KMERS {
            profile[i] = 1.0 / N_KMERS as f64;
        }
    }
    profile
}

/// KL divergence: D_KL(P || Q) = sum P(i) * log(P(i)/Q(i))
fn kl_divergence(p: &[f64; N_KMERS], q: &[f64; N_KMERS]) -> f64 {
    let mut kl = 0.0f64;
    for i in 0..N_KMERS {
        if p[i] > 0.0 && q[i] > 0.0 {
            kl += p[i] * (p[i] / q[i]).ln();
        }
    }
    kl
}

/// GC content of a sequence.
fn gc_content(seq: &[u8]) -> f64 {
    let gc = seq.iter().filter(|&&b| b == b'G' || b == b'C').count();
    let total = seq.iter().filter(|&&b| b == b'A' || b == b'C' || b == b'G' || b == b'T').count();
    if total > 0 { gc as f64 / total as f64 } else { 0.5 }
}

/// Detect prophage regions in a genome.
/// Returns list of (start, end, RegionType) tuples.
pub fn detect_prophage_regions(
    genome: &[u8],
    gene_starts: &[(usize, usize)], // (start, end) of predicted genes
) -> Vec<(usize, usize, RegionType)> {
    let glen = genome.len();
    if glen < WINDOW * 2 {
        return vec![(0, glen, RegionType::Host)];
    }

    // Genome-wide background profiles
    let genome_profile = kmer_profile(genome);
    let genome_gc = gc_content(genome);

    // Pre-compute gene coverage for density calculation
    let mut coverage = vec![false; glen];
    for &(s, e) in gene_starts {
        let start = s.saturating_sub(1).min(glen);
        let end = e.min(glen);
        for i in start..end {
            coverage[i] = true;
        }
    }

    // Sliding window feature extraction
    let mut windows: Vec<WindowFeatures> = Vec::new();
    let mut pos = 0;
    while pos + WINDOW <= glen {
        let win_seq = &genome[pos..pos + WINDOW];

        // k-mer divergence
        let win_profile = kmer_profile(win_seq);
        let kl = kl_divergence(&win_profile, &genome_profile);

        // GC shift
        let local_gc = gc_content(win_seq);
        let gc_shift = (local_gc - genome_gc).abs();

        // Coding density in window
        let covered = coverage[pos..pos + WINDOW].iter().filter(|&&c| c).count();
        let coding_density = covered as f64 / WINDOW as f64;

        // Mean ORF length in window
        let orfs_in_window: Vec<usize> = gene_starts.iter()
            .filter(|&&(s, e)| s >= pos && e <= pos + WINDOW)
            .map(|&(s, e)| e - s)
            .collect();
        let mean_orf_length = if !orfs_in_window.is_empty() {
            orfs_in_window.iter().sum::<usize>() as f64 / orfs_in_window.len() as f64
        } else {
            0.0
        };

        windows.push(WindowFeatures {
            start: pos,
            end: pos + WINDOW,
            kl_divergence: kl,
            coding_density,
            gc_shift,
            mean_orf_length,
            region_type: RegionType::Host, // will be classified below
        });

        pos += STEP;
    }

    if windows.is_empty() {
        return vec![(0, glen, RegionType::Host)];
    }

    // Classify windows using thresholds learned from data
    // Prophage signals: high KL divergence + GC shift + shorter genes + higher density
    let kl_values: Vec<f64> = windows.iter().map(|w| w.kl_divergence).collect();
    let mut sorted_kl = kl_values.clone();
    sorted_kl.sort_by(|a, b| a.partial_cmp(b).unwrap());
    // Prophage threshold: top 5% of KL divergence (most alien windows)
    let kl_threshold = sorted_kl[sorted_kl.len() * 95 / 100];
    // GC shift must be significant
    let gc_threshold = 0.03;

    // Report distribution for calibration
    let kl_median = sorted_kl[sorted_kl.len() / 2];
    let kl_p90 = sorted_kl[sorted_kl.len() * 90 / 100];
    let kl_p95 = sorted_kl[sorted_kl.len() * 95 / 100];
    eprintln!("  KL divergence: median={:.4} p90={:.4} p95={:.4} threshold={:.4}",
        kl_median, kl_p90, kl_p95, kl_threshold);

    for w in windows.iter_mut() {
        // Require BOTH kl divergence AND gc shift for prophage call
        if w.kl_divergence > kl_threshold && w.gc_shift > gc_threshold {
            w.region_type = RegionType::Prophage;
        }
        // Cryptic: very high KL only (top 2%)
        let kl_cryptic = sorted_kl[sorted_kl.len() * 98 / 100];
        if w.kl_divergence > kl_cryptic && w.gc_shift <= gc_threshold {
            w.region_type = RegionType::Cryptic;
        }
    }

    // Merge adjacent windows of same type into regions
    let mut regions: Vec<(usize, usize, RegionType)> = Vec::new();

    // Smooth: require 2+ adjacent prophage windows
    let mut i = 0;
    while i < windows.len() {
        if windows[i].region_type != RegionType::Host {
            // Check if next window is also non-host
            let rt = windows[i].region_type;
            let start = windows[i].start;
            let mut end = windows[i].end;
            let mut count = 1;
            while i + 1 < windows.len() && windows[i + 1].region_type != RegionType::Host {
                i += 1;
                end = windows[i].end;
                count += 1;
            }
            if count >= 2 {
                // Real prophage island (2+ windows = 10kb+)
                regions.push((start, end, rt));
            }
            // else: isolated window, treat as host (noise)
        }
        i += 1;
    }

    if regions.is_empty() {
        return vec![(0, glen, RegionType::Host)];
    }

    // Report
    let total_prophage: usize = regions.iter()
        .filter(|(_, _, t)| *t == RegionType::Prophage)
        .map(|(s, e, _)| e - s).sum();
    let total_cryptic: usize = regions.iter()
        .filter(|(_, _, t)| *t == RegionType::Cryptic)
        .map(|(s, e, _)| e - s).sum();

    for &(rs, re, ref rt) in &regions {
        let label = match rt {
            RegionType::Prophage => "prophage",
            RegionType::Cryptic => "cryptic",
            RegionType::Host => "host",
        };
        eprintln!("  Region: {}..{} ({:.0}kb) {}", rs, re, (re - rs) as f64 / 1000.0, label);
    }
    eprintln!("Prophage detection: {} islands ({:.0}kb prophage, {:.0}kb cryptic)",
        regions.len(), total_prophage as f64 / 1000.0, total_cryptic as f64 / 1000.0);

    regions
}

/// Check if a position falls within a prophage region.
pub fn is_prophage(start: usize, end: usize, regions: &[(usize, usize, RegionType)]) -> bool {
    for &(rs, re, ref rt) in regions {
        if *rt != RegionType::Host && start < re && end > rs {
            return true;
        }
    }
    false
}

/// Check if a position falls within a cryptic (degraded) prophage.
pub fn is_cryptic(start: usize, end: usize, regions: &[(usize, usize, RegionType)]) -> bool {
    for &(rs, re, ref rt) in regions {
        if *rt == RegionType::Cryptic && start < re && end > rs {
            return true;
        }
    }
    false
}
