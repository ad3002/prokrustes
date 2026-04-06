//! Prophage detection via gene-level outlier clustering.
//!
//! Strategy: model the HOST coding well (hex model already trained),
//! then find clusters of consecutive genes with low host-likeness.
//!
//! Key insight from data:
//! - Phage genes have hex_avg ~0.22 vs core ~0.44 (E. coli)
//! - Real phage islands = clusters of ≥5 consecutive low-hex genes
//! - False alarms are typically 3-4 gene clusters (too short)
//! - Island sizes: 10-60kb typically
//!
//! This approach uses predictions from the FIRST pass as input,
//! doesn't need windows or k-mer profiles — just gene-level hex scores.

use crate::types::Gene;

const MIN_CLUSTER_GENES: usize = 5;  // minimum consecutive low-hex genes
const MIN_ISLAND_BP: usize = 8000;   // minimum island size (8kb)
const MAX_GAP_BP: usize = 5000;      // max gap between low-hex genes to keep cluster

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum RegionType {
    Host,
    Prophage,       // cluster of alien genes
    Degraded,       // sparse alien genes (degraded remnant)
}

/// Detect prophage regions from predicted genes using hex outlier clustering.
///
/// Process:
/// 1. Compute genome-specific hex threshold (10th percentile of confident genes)
/// 2. Mark genes below threshold as "outlier"
/// 3. Find clusters of ≥5 consecutive outlier genes
/// 4. Filter by minimum size (8kb)
/// 5. Classify: dense = Prophage, sparse = Degraded
pub fn detect_prophage_regions(
    genes: &[Gene],  // sorted by start position
    genome_len: usize,
) -> Vec<(usize, usize, RegionType)> {
    if genes.len() < 50 {
        return Vec::new();
    }

    // Step 1: Compute adaptive hex threshold from confident genes
    // Use genes ≥600bp as "confident core" — these are almost always real
    let mut confident_hex: Vec<f64> = genes.iter()
        .filter(|g| g.length >= 600 && g.hex_avg > -0.5)
        .map(|g| g.hex_avg)
        .collect();

    if confident_hex.len() < 30 {
        return Vec::new();
    }

    confident_hex.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let threshold = confident_hex[confident_hex.len() / 10]; // 10th percentile

    let median_hex = confident_hex[confident_hex.len() / 2];
    eprintln!("  Prophage detector: host hex median={:.3}, 10th pct={:.3}, threshold={:.3}",
        median_hex, threshold, threshold);

    // Step 2: Mark outlier genes
    let is_outlier: Vec<bool> = genes.iter()
        .map(|g| g.hex_avg < threshold)
        .collect();

    // Step 3: Find clusters of consecutive outlier genes
    // Allow small gaps (1-2 non-outlier genes or physical gap < MAX_GAP_BP)
    let mut clusters: Vec<(usize, usize, usize)> = Vec::new(); // (gene_idx_start, gene_idx_end, n_outliers)
    let mut i = 0;

    while i < genes.len() {
        if !is_outlier[i] {
            i += 1;
            continue;
        }

        // Start a cluster
        let cluster_start = i;
        let mut n_outliers = 1;
        let mut consecutive_non_outlier = 0;
        let mut j = i + 1;

        while j < genes.len() {
            // Check physical gap
            let gap = if genes[j].start > genes[j - 1].end {
                genes[j].start - genes[j - 1].end
            } else {
                0
            };

            if gap > MAX_GAP_BP {
                break; // too far apart, end cluster
            }

            if is_outlier[j] {
                n_outliers += 1;
                consecutive_non_outlier = 0;
            } else {
                consecutive_non_outlier += 1;
                if consecutive_non_outlier >= 3 {
                    // 3+ non-outlier genes in a row = cluster break
                    j -= consecutive_non_outlier - 1;
                    break;
                }
            }
            j += 1;
        }

        let cluster_end = j - 1;
        if n_outliers >= MIN_CLUSTER_GENES {
            clusters.push((cluster_start, cluster_end, n_outliers));
        }

        i = j;
    }

    // Step 4: Convert gene clusters to genomic regions, filter by size
    let mut regions: Vec<(usize, usize, RegionType)> = Vec::new();

    for (gs, ge, n_outliers) in &clusters {
        let start = genes[*gs].start.saturating_sub(1000); // extend 1kb
        let end = (genes[*ge].end + 1000).min(genome_len);
        let size = end - start;

        if size < MIN_ISLAND_BP {
            continue;
        }

        let n_genes = ge - gs + 1;
        let outlier_frac = *n_outliers as f64 / n_genes as f64;
        let gene_density = n_genes as f64 / (size as f64 / 1000.0);

        // Classify
        let rtype = if outlier_frac > 0.5 && gene_density > 1.0 {
            RegionType::Prophage  // dense cluster of alien genes
        } else {
            RegionType::Degraded  // sparse, likely degraded remnant
        };

        let label = match rtype {
            RegionType::Prophage => "prophage",
            RegionType::Degraded => "degraded",
            RegionType::Host => "host",
        };
        eprintln!("  Island: {}..{} ({:.0}kb, {} genes, {} outliers, {:.1} genes/kb) {}",
            start, end, size as f64 / 1000.0, n_genes, n_outliers, gene_density, label);

        regions.push((start, end, rtype));
    }

    if !regions.is_empty() {
        let n_prophage = regions.iter().filter(|(_, _, t)| *t == RegionType::Prophage).count();
        let kb_prophage: usize = regions.iter()
            .filter(|(_, _, t)| *t == RegionType::Prophage)
            .map(|(s, e, _)| e - s).sum();
        let n_degraded = regions.iter().filter(|(_, _, t)| *t == RegionType::Degraded).count();
        eprintln!("Prophage detection: {} prophage ({:.0}kb), {} degraded",
            n_prophage, kb_prophage as f64 / 1000.0, n_degraded);
    }

    regions
}

/// Compute region bonus for a gene based on overlap with detected regions.
pub fn region_bonus_for(
    start: usize,
    end: usize,
    regions: &[(usize, usize, RegionType)],
) -> f64 {
    for &(rs, re, ref rt) in regions {
        if start < re && end > rs {
            let ov = end.min(re) - start.max(rs);
            let gene_len = end - start;
            if ov * 2 >= gene_len {
                return match rt {
                    RegionType::Prophage => 0.0,   // no bonus yet — validate first
                    RegionType::Degraded => -0.04,  // penalize degraded (pseudogene-rich)
                    RegionType::Host => 0.0,
                };
            }
        }
    }
    0.0
}
