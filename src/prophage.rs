//! Prophage detection via gene-level scoring law.
//!
//! Score(S) = w1·lowHexFrac + w2·RunScore + w3·densityExcess - w4·fragmentation
//!
//! where RunScore = a1·maxRun/core99 + a2·runMass
//!
//! No penalty for ameliorated islands. Only penalty for incoherent signal.

use crate::types::Gene;

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum RegionType {
    Host,
    Prophage,
    Degraded,
}

/// Summary features for a candidate segment.
#[derive(Debug)]
struct SegmentFeatures {
    start: usize,
    end: usize,
    n_genes: usize,
    span_kb: f64,
    low_hex_frac: f64,
    max_run: usize,
    run_mass: f64,       // fraction of low-hex genes in runs ≥3
    fragmentation: f64,  // n_runs / n_low_hex_genes (1.0 = all singletons)
    density: f64,        // genes/kb in segment
    density_excess: f64, // (density - core_density) / core_density
}

/// Compute genome-specific background statistics.
struct GenomeStats {
    hex_threshold: f64,  // 10th percentile of confident genes
    core_density: f64,   // genes/kb in core regions
    core99_run: usize,   // 99th percentile of low-hex run length in core
    median_hex: f64,
}

fn compute_genome_stats(genes: &[Gene], genome_len: usize) -> Option<GenomeStats> {
    if genes.len() < 50 { return None; }

    // Hex threshold: 10th percentile of confident (≥600bp) genes
    let mut confident_hex: Vec<f64> = genes.iter()
        .filter(|g| g.length >= 600 && g.hex_avg > -1.0)
        .map(|g| g.hex_avg)
        .collect();
    if confident_hex.len() < 30 { return None; }
    confident_hex.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let threshold = confident_hex[confident_hex.len() / 10];
    let median_hex = confident_hex[confident_hex.len() / 2];

    // Core density: total genes / genome size
    let core_density = genes.len() as f64 / (genome_len as f64 / 1000.0);

    // Core99 run: compute all low-hex runs, take 99th percentile
    let runs = find_low_hex_runs(genes, threshold);
    let mut run_lengths: Vec<usize> = runs.iter().map(|r| r.len).collect();
    run_lengths.sort();
    let core99 = if run_lengths.is_empty() {
        3
    } else {
        run_lengths[run_lengths.len() * 99 / 100].max(3)
    };

    eprintln!("  Prophage stats: hex_threshold={:.3}, core_density={:.2} g/kb, core99_run={}",
        threshold, core_density, core99);

    Some(GenomeStats { hex_threshold: threshold, core_density, core99_run: core99, median_hex })
}

struct LowHexRun {
    start_idx: usize, // index in genes array
    end_idx: usize,   // inclusive
    len: usize,
}

/// Find all runs of consecutive low-hex genes (allowing 1-2 gene gaps).
fn find_low_hex_runs(genes: &[Gene], threshold: f64) -> Vec<LowHexRun> {
    let is_low: Vec<bool> = genes.iter().map(|g| g.hex_avg < threshold).collect();
    let mut runs = Vec::new();
    let mut i = 0;

    while i < genes.len() {
        if !is_low[i] { i += 1; continue; }

        let start = i;
        let mut low_count = 1;
        let mut consecutive_high = 0;
        let mut j = i + 1;

        while j < genes.len() {
            // Check physical gap
            let gap = if genes[j].start > genes[j-1].end {
                genes[j].start - genes[j-1].end
            } else { 0 };
            if gap > 8000 { break; } // >8kb gap = cluster break

            if is_low[j] {
                low_count += 1;
                consecutive_high = 0;
            } else {
                consecutive_high += 1;
                if consecutive_high >= 3 { break; } // 3+ high in row = break
            }
            j += 1;
        }

        // Trim trailing high genes
        while j > start + 1 && !is_low[j - 1] { j -= 1; }

        runs.push(LowHexRun { start_idx: start, end_idx: j - 1, len: low_count });
        i = j;
    }

    runs
}

/// Compute segment features for a candidate region (defined by gene range).
fn compute_segment_features(
    genes: &[Gene],
    start_idx: usize,
    end_idx: usize,
    stats: &GenomeStats,
    all_runs: &[LowHexRun],
) -> SegmentFeatures {
    let seg_genes = &genes[start_idx..=end_idx];
    let n_genes = seg_genes.len();
    let start = seg_genes[0].start;
    let end = seg_genes[seg_genes.len() - 1].end;
    let span_kb = (end - start) as f64 / 1000.0;

    // lowHexFrac
    let n_low = seg_genes.iter().filter(|g| g.hex_avg < stats.hex_threshold).count();
    let low_hex_frac = n_low as f64 / n_genes.max(1) as f64;

    // Run analysis within this segment
    let seg_runs: Vec<&LowHexRun> = all_runs.iter()
        .filter(|r| r.start_idx >= start_idx && r.end_idx <= end_idx)
        .collect();

    let max_run = seg_runs.iter().map(|r| r.len).max().unwrap_or(0);

    // runMass: fraction of low-hex genes that are in runs ≥3
    let genes_in_runs3: usize = seg_runs.iter()
        .filter(|r| r.len >= 3)
        .map(|r| r.len)
        .sum();
    let run_mass = if n_low > 0 { genes_in_runs3 as f64 / n_low as f64 } else { 0.0 };

    // fragmentation: n_runs / n_low_hex (1.0 = all singletons, 0.05 = one big block)
    let n_runs = seg_runs.len();
    let fragmentation = if n_low > 0 { n_runs as f64 / n_low as f64 } else { 1.0 };

    // density
    let density = if span_kb > 0.0 { n_genes as f64 / span_kb } else { 0.0 };
    let density_excess = if stats.core_density > 0.0 {
        ((density - stats.core_density) / stats.core_density).max(0.0)
    } else { 0.0 };

    SegmentFeatures {
        start, end, n_genes, span_kb, low_hex_frac,
        max_run, run_mass, fragmentation, density, density_excess,
    }
}

/// Score a segment using the scoring law.
fn score_segment(feat: &SegmentFeatures, stats: &GenomeStats) -> f64 {
    let w1 = 1.0;  // lowHexFrac weight
    let w2 = 0.6;  // RunScore weight
    let w3 = 0.3;  // densityExcess weight
    let w4 = 0.5;  // fragmentation penalty weight

    let a1 = 0.7;  // maxRun component of RunScore
    let a2 = 0.3;  // runMass component of RunScore

    // F(S) = clip(lowHexFrac / 0.5, 0, 1)  — saturates at 50%
    let f = (feat.low_hex_frac / 0.5).min(1.0).max(0.0);

    // RunScore = a1 * maxRun/core99 + a2 * runMass, clipped to [0, 2]
    let core99 = stats.core99_run.max(1) as f64;
    let run_score = (a1 * feat.max_run as f64 / core99 + a2 * feat.run_mass).min(2.0).max(0.0);

    // D(S) = positive density excess only, clipped to [0, 1]
    let d = feat.density_excess.min(1.0);

    // C(S) = fragmentation penalty: high if signal is scattered singletons
    let c = if feat.fragmentation > 0.8 { 1.0 }
            else if feat.fragmentation > 0.5 { 0.5 }
            else { 0.0 };

    w1 * f + w2 * run_score + w3 * d - w4 * c
}

/// Main entry point: detect prophage regions.
pub fn detect_prophage_regions(
    genes: &[Gene],
    genome_len: usize,
) -> Vec<(usize, usize, RegionType)> {
    let stats = match compute_genome_stats(genes, genome_len) {
        Some(s) => s,
        None => return Vec::new(),
    };

    let all_runs = find_low_hex_runs(genes, stats.hex_threshold);

    // Generate candidate segments from runs
    // Merge nearby runs (within 10 genes) into candidate segments
    let mut candidates: Vec<(usize, usize)> = Vec::new(); // (start_idx, end_idx)

    let significant_runs: Vec<&LowHexRun> = all_runs.iter()
        .filter(|r| r.len >= 3)
        .collect();

    if significant_runs.is_empty() {
        return Vec::new();
    }

    // Merge overlapping/adjacent significant runs
    let mut i = 0;
    while i < significant_runs.len() {
        let mut seg_start = significant_runs[i].start_idx;
        let mut seg_end = significant_runs[i].end_idx;

        // Extend: merge with next run if gap < 10 genes
        while i + 1 < significant_runs.len() {
            let next = significant_runs[i + 1];
            if next.start_idx <= seg_end + 10 {
                seg_end = seg_end.max(next.end_idx);
                i += 1;
            } else {
                break;
            }
        }

        // Extend segment boundaries by 2 genes on each side
        seg_start = seg_start.saturating_sub(2);
        seg_end = (seg_end + 2).min(genes.len() - 1);

        candidates.push((seg_start, seg_end));
        i += 1;
    }

    // Score each candidate
    let mut regions: Vec<(usize, usize, RegionType)> = Vec::new();
    let score_threshold = 0.6; // minimum score to call foreign

    for (seg_start, seg_end) in &candidates {
        let feat = compute_segment_features(genes, *seg_start, *seg_end, &stats, &all_runs);

        // Minimum size filter
        if feat.span_kb < 5.0 || feat.n_genes < 4 { continue; }

        let score = score_segment(&feat, &stats);

        if score >= score_threshold {
            let rtype = if feat.low_hex_frac > 0.35 && feat.max_run >= 3 {
                RegionType::Prophage
            } else {
                RegionType::Degraded
            };

            let label = match rtype {
                RegionType::Prophage => "prophage",
                RegionType::Degraded => "degraded",
                RegionType::Host => "host",
            };

            eprintln!("  Island: {}..{} ({:.0}kb, {} genes) score={:.2} [F={:.2} R={:.2} D={:.2} frag={:.2}] {}",
                feat.start, feat.end, feat.span_kb, feat.n_genes,
                score, feat.low_hex_frac, feat.run_mass, feat.density_excess, feat.fragmentation,
                label);

            regions.push((feat.start, feat.end, rtype));
        }
    }

    if !regions.is_empty() {
        let n_prophage = regions.iter().filter(|(_, _, t)| *t == RegionType::Prophage).count();
        let kb_total: usize = regions.iter().map(|(s, e, _)| e - s).sum();
        eprintln!("Prophage detection: {} regions ({:.0}kb total, {} prophage)",
            regions.len(), kb_total as f64 / 1000.0, n_prophage);
    }

    regions
}

/// Compute region bonus for routing.
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
                    RegionType::Prophage => 0.0,
                    RegionType::Degraded => -0.03,
                    RegionType::Host => 0.0,
                };
            }
        }
    }
    0.0
}
