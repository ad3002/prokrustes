//! Multi-channel prophage detection and region routing.
//!
//! Architecture: 4-layer cascade
//! Layer 1: 5-channel window features (composition + ORF architecture)
//! Layer 2: One-class anomaly scoring (robust Z-scores, multi-channel convergence)
//! Layer 3: Region segmentation (smoothing, merging, classification)
//! Layer 4: Routing via threshold adjustment (region_bonus on Gene)
//!
//! Key principle: model the HOST well (one-class), flag anything that deviates
//! across multiple channels simultaneously. This naturally handles prophage,
//! IS-rich regions, and HGT islands without phage-specific signatures.

use crate::io::nt4;

const WINDOW: usize = 5000;   // 5kb sliding window
const STEP: usize = 500;      // 500bp step
const KMER_K: usize = 4;      // tetranucleotide
const N_KMERS: usize = 256;   // 4^4
const MIN_REGION_KB: usize = 10000; // minimum island size 10kb
const SIGMA_THRESHOLD: f64 = 2.5;   // anomaly threshold in sigma units
const MIN_CHANNELS: usize = 2;      // minimum anomalous channels for foreignness

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum RegionType {
    Host,
    Prophage,         // structured prophage (alien composition + compact architecture)
    Degraded,         // degraded prophage (alien composition + low coding density)
    MobileNonphage,   // IS-rich island, HGT region
}

/// Per-window feature vector.
pub struct WindowFeatures {
    pub start: usize,
    pub end: usize,
    // 5 channels
    pub tetra_distance: f64,   // tetranucleotide distance from host (sigma units)
    pub gc_deviation: f64,     // |local_gc - genome_gc| (sigma units)
    pub coding_density: f64,   // fraction of window covered by ORFs
    pub short_orf_frac: f64,   // fraction of ORFs <300bp in window
    pub overlap_frac: f64,     // fraction of window in overlapping ORFs
    // Derived
    pub foreignness: f64,      // composite anomaly score
    pub region_type: RegionType,
}

type TetraProfile = [f64; N_KMERS];

// ---- Layer 1: Feature computation ----

/// Compute tetranucleotide frequency profile for a sequence.
fn tetra_profile(seq: &[u8]) -> TetraProfile {
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
        if valid { counts[idx] += 1; total += 1; }
    }
    let mut profile = [0.0f64; N_KMERS];
    let t = total as f64;
    if t > 0.0 {
        for i in 0..N_KMERS {
            profile[i] = counts[i] as f64 / t;
        }
    }
    profile
}

/// Euclidean distance between two tetra profiles.
fn tetra_distance(a: &TetraProfile, b: &TetraProfile) -> f64 {
    let mut sum = 0.0f64;
    for i in 0..N_KMERS {
        let d = a[i] - b[i];
        sum += d * d;
    }
    sum.sqrt()
}

/// GC content of a sequence.
fn gc_content(seq: &[u8]) -> f64 {
    let gc = seq.iter().filter(|&&b| b == b'G' || b == b'C').count();
    let total = seq.iter().filter(|&&b| matches!(b, b'A' | b'C' | b'G' | b'T')).count();
    if total > 0 { gc as f64 / total as f64 } else { 0.5 }
}

/// Compute all 5-channel features for sliding windows.
/// gene_coords: (start, end, length) from first-pass predictions.
pub fn compute_window_features(
    genome: &[u8],
    gene_coords: &[(usize, usize, usize)], // (start, end, length)
) -> Vec<WindowFeatures> {
    let glen = genome.len();
    if glen < WINDOW * 2 { return Vec::new(); }

    // Genome-wide references
    let genome_tetra = tetra_profile(genome);
    let genome_gc = gc_content(genome);

    // Pre-compute per-base coverage and overlap arrays
    let mut coverage = vec![0u8; glen]; // 0=uncovered, 1=covered, 2+=overlap
    for &(s, e, _) in gene_coords {
        let start = s.saturating_sub(1).min(glen);
        let end = e.min(glen);
        for i in start..end {
            coverage[i] = coverage[i].saturating_add(1);
        }
    }

    // Sliding window feature extraction
    let mut windows: Vec<WindowFeatures> = Vec::new();
    let mut pos = 0;
    while pos + WINDOW <= glen {
        let win_seq = &genome[pos..pos + WINDOW];

        // Channel 1: tetranucleotide distance
        let win_tetra = tetra_profile(win_seq);
        let td = tetra_distance(&win_tetra, &genome_tetra);

        // Channel 2: GC deviation
        let local_gc = gc_content(win_seq);
        let gc_dev = (local_gc - genome_gc).abs();

        // Channel 3: coding density
        let covered = coverage[pos..pos + WINDOW].iter().filter(|&&c| c >= 1).count();
        let coding_dens = covered as f64 / WINDOW as f64;

        // Channel 4: short ORF fraction
        let orfs_in_window: Vec<usize> = gene_coords.iter()
            .filter(|&&(s, e, _)| s >= pos && e <= pos + WINDOW)
            .map(|&(_, _, l)| l)
            .collect();
        let n_orfs = orfs_in_window.len();
        let n_short = orfs_in_window.iter().filter(|&&l| l < 300).count();
        let short_frac = if n_orfs > 0 { n_short as f64 / n_orfs as f64 } else { 0.0 };

        // Channel 5: overlap fraction
        let overlapping = coverage[pos..pos + WINDOW].iter().filter(|&&c| c >= 2).count();
        let overlap_frac = overlapping as f64 / WINDOW as f64;

        windows.push(WindowFeatures {
            start: pos,
            end: pos + WINDOW,
            tetra_distance: td,
            gc_deviation: gc_dev,
            coding_density: coding_dens,
            short_orf_frac: short_frac,
            overlap_frac: overlap_frac,
            foreignness: 0.0,
            region_type: RegionType::Host,
        });

        pos += STEP;
    }

    windows
}

// ---- Layer 2: One-class anomaly scoring ----

/// Robust Z-scores using median and MAD (median absolute deviation).
/// Returns Z-scores in sigma units. Resistant to outliers.
fn robust_zscores(values: &[f64]) -> Vec<f64> {
    if values.is_empty() { return Vec::new(); }

    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = sorted[sorted.len() / 2];

    let mut abs_devs: Vec<f64> = values.iter().map(|&v| (v - median).abs()).collect();
    abs_devs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mad = abs_devs[abs_devs.len() / 2];
    let sigma = 1.4826 * mad; // MAD to σ conversion for normal distribution

    if sigma < 1e-10 {
        return vec![0.0; values.len()];
    }

    values.iter().map(|&v| (v - median).abs() / sigma).collect()
}

/// Compute foreignness: multi-channel anomaly score.
/// A window is foreign when ≥MIN_CHANNELS channels exceed SIGMA_THRESHOLD.
pub fn compute_foreignness(windows: &mut [WindowFeatures]) {
    if windows.is_empty() { return; }

    // Extract each channel
    let tetra_vals: Vec<f64> = windows.iter().map(|w| w.tetra_distance).collect();
    let gc_vals: Vec<f64> = windows.iter().map(|w| w.gc_deviation).collect();
    // For coding density and short_orf: anomaly = deviation from TYPICAL, not just high
    let cd_vals: Vec<f64> = windows.iter().map(|w| w.coding_density).collect();
    let so_vals: Vec<f64> = windows.iter().map(|w| w.short_orf_frac).collect();
    let ov_vals: Vec<f64> = windows.iter().map(|w| w.overlap_frac).collect();

    // Robust Z-scores per channel
    let z_tetra = robust_zscores(&tetra_vals);
    let z_gc = robust_zscores(&gc_vals);
    let z_cd = robust_zscores(&cd_vals);
    let z_so = robust_zscores(&so_vals);
    let z_ov = robust_zscores(&ov_vals);

    for (i, w) in windows.iter_mut().enumerate() {
        // Count channels exceeding threshold
        let mut n_anomalous = 0u32;
        let mut excess_sum = 0.0f64;

        let channels = [z_tetra[i], z_gc[i], z_cd[i], z_so[i], z_ov[i]];
        for &z in &channels {
            if z > SIGMA_THRESHOLD {
                n_anomalous += 1;
                excess_sum += z - SIGMA_THRESHOLD;
            }
        }

        // Multi-channel convergence: require ≥2 anomalous channels
        if n_anomalous >= MIN_CHANNELS as u32 {
            w.foreignness = excess_sum;
        }
    }
}

// ---- Layer 3: Segmentation ----

/// Segment windows into regions.
pub fn segment_regions(
    windows: &[WindowFeatures],
    genome_len: usize,
) -> Vec<(usize, usize, RegionType)> {
    if windows.is_empty() {
        return vec![];
    }

    // Step 1: 3-window median filter on foreignness
    let mut smoothed: Vec<f64> = vec![0.0; windows.len()];
    for i in 0..windows.len() {
        let mut vals = Vec::new();
        if i > 0 { vals.push(windows[i - 1].foreignness); }
        vals.push(windows[i].foreignness);
        if i + 1 < windows.len() { vals.push(windows[i + 1].foreignness); }
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        smoothed[i] = vals[vals.len() / 2];
    }

    // Step 2: Find runs of anomalous windows
    let mut regions: Vec<(usize, usize, RegionType)> = Vec::new();
    let mut i = 0;
    while i < windows.len() {
        if smoothed[i] > 0.0 {
            let start = windows[i].start;
            let mut end = windows[i].end;
            let mut sum_foreign = smoothed[i];
            let mut count = 1;

            // Extend: allow 1 gap window (bridge small holes)
            while i + 1 < windows.len() {
                if smoothed[i + 1] > 0.0 {
                    i += 1;
                    end = windows[i].end;
                    sum_foreign += smoothed[i];
                    count += 1;
                } else if i + 2 < windows.len() && smoothed[i + 2] > 0.0 {
                    // Bridge 1 gap
                    i += 2;
                    end = windows[i].end;
                    sum_foreign += smoothed[i];
                    count += 1;
                } else {
                    break;
                }
            }

            let size = end - start;
            if size >= MIN_REGION_KB {
                // Classify based on feature profile within region
                let region_windows: Vec<&WindowFeatures> = windows.iter()
                    .filter(|w| w.start >= start && w.end <= end)
                    .collect();

                let avg_cd: f64 = region_windows.iter().map(|w| w.coding_density).sum::<f64>()
                    / region_windows.len().max(1) as f64;
                let avg_so: f64 = region_windows.iter().map(|w| w.short_orf_frac).sum::<f64>()
                    / region_windows.len().max(1) as f64;
                let avg_tetra: f64 = region_windows.iter().map(|w| w.tetra_distance).sum::<f64>()
                    / region_windows.len().max(1) as f64;

                let rtype = if avg_cd < 0.65 {
                    RegionType::Degraded // low coding = degraded remnant
                } else if avg_so > 0.35 || avg_cd > 0.90 {
                    RegionType::Prophage // compact + many short genes = structured phage
                } else {
                    RegionType::MobileNonphage // alien but normal architecture
                };

                regions.push((start, end, rtype));
            }
        }
        i += 1;
    }

    // Report
    if !regions.is_empty() {
        let mut n_prophage = 0usize;
        let mut n_degraded = 0usize;
        let mut n_mobile = 0usize;
        let mut kb_prophage = 0usize;
        let mut kb_degraded = 0usize;

        for &(s, e, ref rt) in &regions {
            let label = match rt {
                RegionType::Prophage => { n_prophage += 1; kb_prophage += e - s; "prophage" },
                RegionType::Degraded => { n_degraded += 1; kb_degraded += e - s; "degraded" },
                RegionType::MobileNonphage => { n_mobile += 1; "mobile" },
                RegionType::Host => "host",
            };
            eprintln!("  Region: {}..{} ({:.0}kb) {}", s, e, (e - s) as f64 / 1000.0, label);
        }
        eprintln!("Prophage detection: {} prophage ({:.0}kb), {} degraded ({:.0}kb), {} mobile",
            n_prophage, kb_prophage as f64 / 1000.0,
            n_degraded, kb_degraded as f64 / 1000.0,
            n_mobile);
    }

    regions
}

// ---- Layer 4: Routing ----

/// Compute region bonus for a gene based on its overlap with detected regions.
pub fn region_bonus_for(
    start: usize,
    end: usize,
    regions: &[(usize, usize, RegionType)],
) -> f64 {
    for &(rs, re, ref rt) in regions {
        // Check overlap
        if start < re && end > rs {
            let ov = end.min(re) - start.max(rs);
            let gene_len = end - start;
            // Only apply if >50% of gene is in the region
            if ov * 2 >= gene_len {
                return match rt {
                    RegionType::Prophage => 0.0,         // no bonus (avoid adding FP)
                    RegionType::Degraded => -0.05,       // penalize degraded (pseudogene-rich)
                    RegionType::MobileNonphage => 0.0,
                    RegionType::Host => 0.0,
                };
            }
        }
    }
    0.0
}

// ---- Public entry point ----

/// Full prophage detection pipeline: features → anomaly scoring → segmentation.
pub fn detect_prophage_regions(
    genome: &[u8],
    gene_coords: &[(usize, usize, usize)],
) -> Vec<(usize, usize, RegionType)> {
    let mut windows = compute_window_features(genome, gene_coords);
    if windows.is_empty() {
        return Vec::new();
    }
    compute_foreignness(&mut windows);
    segment_regions(&windows, genome.len())
}
