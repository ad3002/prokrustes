use std::collections::{HashMap, HashSet};

use crate::types::Gene;

pub fn dp_select(genes: &[Gene], max_ov: usize) -> Vec<usize> {
    if genes.is_empty() { return vec![]; }

    let n = genes.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| genes[a].end.cmp(&genes[b].end)
        .then_with(|| genes[b].weight.partial_cmp(&genes[a].weight)
            .unwrap_or(std::cmp::Ordering::Equal)));

    let ends: Vec<usize> = idx.iter().map(|&i| genes[i].end).collect();

    let mut p = vec![-1i64; n];
    for i in 0..n {
        let gi = idx[i];
        let target = genes[gi].start + max_ov - 1;
        let mut lo: isize = 0;
        let mut hi: isize = i as isize - 1;
        let mut result: i64 = -1;
        while lo <= hi {
            let mid = (lo + hi) as usize / 2;
            if ends[mid] <= target {
                result = mid as i64;
                lo = mid as isize + 1;
            } else {
                hi = mid as isize - 1;
            }
        }
        p[i] = result;
    }

    // Forward DP
    let mut dp = vec![0.0f64; n + 1];
    for i in 1..=n {
        let gi = idx[i-1];
        let mut w = genes[gi].weight;
        let sp = genes[gi].shadow_pen;
        if sp < 1.0 { w *= sp; }
        let prev = (p[i-1] + 1) as usize;
        dp[i] = dp[i-1].max(dp[prev] + w);
    }

    // Traceback
    let mut sel = Vec::new();
    let mut i = n;
    while i > 0 {
        let gi = idx[i-1];
        let mut w = genes[gi].weight;
        let sp = genes[gi].shadow_pen;
        if sp < 1.0 { w *= sp; }
        let prev = (p[i-1] + 1) as usize;
        if (dp[i] - (dp[prev] + w)).abs() < 1e-10 {
            sel.push(gi);
            i = prev;
        } else {
            i -= 1;
        }
    }
    sel
}

pub fn iterative_dp(genes: &mut [Gene], max_ov: usize, n_iter: usize) -> Vec<usize> {
    for g in genes.iter_mut() {
        g._base_weight = g.weight;
        g.adj_bonus = 0.0;
    }
    let mut result = Vec::new();

    for it in 0..n_iter {
        for g in genes.iter_mut() {
            g.weight = g._base_weight + g.adj_bonus;
        }
        result = dp_select(genes, max_ov);

        if it < n_iter - 1 {
            let mut sel_sorted: Vec<usize> = result.clone();
            sel_sorted.sort_by_key(|&i| genes[i].start);
            let sel_starts: Vec<usize> = sel_sorted.iter().map(|&i| genes[i].start).collect();

            for ci in 0..genes.len() {
                let mut bonus = 0.0f64;
                let target_lo = genes[ci].start.saturating_sub(300);
                let idx = sel_starts.partition_point(|&s| s < target_lo);

                for k in idx..sel_sorted.len().min(idx + 30) {
                    let si = sel_sorted[k];
                    if genes[si].start > genes[ci].end + 300 { break; }
                    if genes[si].is_plus != genes[ci].is_plus { continue; }
                    if genes[si].start == genes[ci].start && genes[si].end == genes[ci].end { continue; }

                    let gap: isize;
                    if genes[si].end <= genes[ci].start {
                        gap = (genes[ci].start - genes[si].end) as isize;
                    } else if genes[ci].end <= genes[si].start {
                        gap = (genes[si].start - genes[ci].end) as isize;
                    } else {
                        let ov = genes[si].end.min(genes[ci].end) as isize
                            - genes[si].start.max(genes[ci].start) as isize + 1;
                        if ov <= 4 { gap = -ov; } else { continue; }
                    }

                    let b = if gap >= -4 && gap <= 4 { 0.15 }
                        else if gap <= 30 { 0.10 }
                        else if gap <= 80 { 0.05 }
                        else if gap <= 150 { 0.02 }
                        else { 0.0 };
                    bonus = bonus.max(b);
                }
                genes[ci].adj_bonus = bonus;
            }
        }
    }
    result
}

pub fn gap_fill(selected: &mut Vec<usize>, all: &[Gene], candidates: &[usize], min_score: f64, min_len: usize) {
    let sel_set: HashSet<(usize, usize, bool)> = selected.iter().map(|&i| {
        (all[i].start, all[i].end, all[i].is_plus)
    }).collect();

    let mut remaining: Vec<usize> = candidates.iter().filter(|&&i| {
        !sel_set.contains(&(all[i].start, all[i].end, all[i].is_plus))
    }).copied().collect();
    remaining.sort_by(|&a, &b| all[b].score.partial_cmp(&all[a].score).unwrap_or(std::cmp::Ordering::Equal));

    for &oi in &remaining {
        let o = &all[oi];
        if o.score < min_score || o.length < min_len { continue; }
        let o_len = (o.end - o.start + 1) as f64;
        let mut ok = true;
        for &ri in selected.iter() {
            let a = &all[ri];
            let ov_s = o.start.max(a.start);
            let ov_e = o.end.min(a.end);
            if ov_e >= ov_s {
                let overlap = (ov_e - ov_s + 1) as f64;
                let max_ov = if o.is_plus == a.is_plus { 45.0 } else { 90.0 };
                if overlap > max_ov || overlap / o_len > 0.25 {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            selected.push(oi);
        }
    }
}

/// Targeted gap rescue (cycle 20): find large intergenic gaps and try to fill them.
/// Looks for gaps > min_gap bp in selected genes and inserts candidates that fit.
pub fn gap_targeted_rescue(selected: &mut Vec<usize>, all: &[Gene], genome_len: usize, min_gap: usize) {
    // Collect selected gene positions
    let mut positions: Vec<(usize, usize)> = selected.iter()
        .map(|&i| (all[i].start, all[i].end))
        .collect();
    positions.sort();

    // Find large gaps
    let mut gaps: Vec<(usize, usize)> = Vec::new();
    let mut prev_end = 0;
    for &(start, end) in &positions {
        if start > prev_end + min_gap {
            gaps.push((prev_end, start));
        }
        prev_end = prev_end.max(end);
    }
    if genome_len > prev_end + min_gap {
        gaps.push((prev_end, genome_len));
    }

    if gaps.is_empty() { return; }

    let sel_set: HashSet<(usize, usize, bool)> = selected.iter()
        .map(|&i| (all[i].start, all[i].end, all[i].is_plus)).collect();

    let mut added = 0;
    for (gap_start, gap_end) in &gaps {
        let mut gap_cands: Vec<(usize, f64)> = Vec::new();
        for (i, o) in all.iter().enumerate() {
            if sel_set.contains(&(o.start, o.end, o.is_plus)) { continue; }
            if o.start >= *gap_start && o.end <= *gap_end
                && o.score > 0.42 && o.length >= 150 && o.is_longest && o.is_atg() {
                gap_cands.push((i, o.score));
            }
        }
        gap_cands.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

        for &(idx, _) in gap_cands.iter().take(5) {
            let o = &all[idx];
            let o_len = (o.end - o.start + 1) as f64;
            let mut ok = true;
            for &si in selected.iter() {
                let a = &all[si];
                let ov_s = o.start.max(a.start);
                let ov_e = o.end.min(a.end);
                if ov_e >= ov_s {
                    let overlap = (ov_e - ov_s + 1) as f64;
                    if overlap > 30.0 || overlap / o_len > 0.15 {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                selected.push(idx);
                added += 1;
            }
        }
    }
}

/// Filter same-strand overlaps > max_ov bp (cycle 20, post-DP cleanup).
/// If two genes on same strand overlap too much, remove the lower-scoring one.
pub fn filter_same_strand_overlaps(selected: &mut Vec<usize>, all: &[Gene], max_ov: usize) {
    let mut sorted: Vec<usize> = selected.clone();
    sorted.sort_by_key(|&i| all[i].start);
    let mut keep = vec![true; sorted.len()];
    for i in 0..sorted.len() {
        if !keep[i] { continue; }
        let a = &all[sorted[i]];
        for j in (i + 1)..sorted.len() {
            if !keep[j] { continue; }
            let b = &all[sorted[j]];
            if b.start > a.end { break; }
            if a.is_plus != b.is_plus { continue; }
            let ov_s = a.start.max(b.start);
            let ov_e = a.end.min(b.end);
            if ov_e >= ov_s {
                let overlap = ov_e - ov_s + 1;
                if overlap > max_ov {
                    // Remove lower-scoring
                    if a.score >= b.score { keep[j] = false; }
                    else { keep[i] = false; break; }
                }
            }
        }
    }
    *selected = sorted.into_iter().enumerate()
        .filter(|&(i, _)| keep[i]).map(|(_, v)| v).collect();
}

/// Filter opposite-strand overlaps > max_ov bp.
/// Prodigal: max 200bp opposite-strand overlap, and < 50% of either gene.
/// Remove the shorter gene in conflicting pairs.
pub fn filter_opposite_strand_overlaps(selected: &mut Vec<usize>, all: &[Gene], max_ov: usize) {
    let mut sorted: Vec<usize> = selected.clone();
    sorted.sort_by_key(|&i| all[i].start);
    let mut keep = vec![true; sorted.len()];
    for i in 0..sorted.len() {
        if !keep[i] { continue; }
        let a = &all[sorted[i]];
        for j in (i + 1)..sorted.len() {
            if !keep[j] { continue; }
            let b = &all[sorted[j]];
            if b.start > a.end { break; }
            // Only opposite strand
            if a.is_plus == b.is_plus { continue; }
            let ov_s = a.start.max(b.start);
            let ov_e = a.end.min(b.end);
            if ov_e >= ov_s {
                let overlap = ov_e - ov_s + 1;
                let a_len = a.end - a.start + 1;
                let b_len = b.end - b.start + 1;
                // Reject if overlap > max_ov OR > 50% of either gene
                if overlap > max_ov || overlap * 2 > a_len || overlap * 2 > b_len {
                    // Remove the shorter/lower-scoring gene
                    if a.length > b.length || (a.length == b.length && a.score >= b.score) {
                        keep[j] = false;
                    } else {
                        keep[i] = false;
                        break;
                    }
                }
            }
        }
    }
    let removed = keep.iter().filter(|&&k| !k).count();
    if removed > 0 {
        eprintln!("Opposite-strand overlap filter: {} genes removed", removed);
    }
    *selected = sorted.into_iter().enumerate()
        .filter(|&(i, _)| keep[i]).map(|(_, v)| v).collect();
}

pub fn detect_shadows(genes: &mut [Gene]) {
    let n = genes.len();
    if n < 2 { return; }
    let mut sorted: Vec<usize> = (0..n).collect();
    sorted.sort_by_key(|&i| genes[i].start);

    for ii in 0..n {
        let i = sorted[ii];
        let gs = genes[i].start;
        let ge = genes[i].end;
        let gstrand = genes[i].is_plus;
        let glen = genes[i].length;
        let gcod = genes[i].hex_avg;

        let lo = if ii >= 50 { ii - 50 } else { 0 };
        let hi = (ii + 50).min(n);

        for jj in lo..hi {
            if jj == ii { continue; }
            let j = sorted[jj];
            if genes[j].is_plus == gstrand { continue; }

            let ov_s = gs.max(genes[j].start);
            let ov_e = ge.min(genes[j].end);
            if ov_e < ov_s { continue; }

            let overlap = (ov_e - ov_s + 1) as f64;
            let orf_len = (ge - gs + 1) as f64;
            let other_cod = genes[j].hex_avg;
            let other_len = genes[j].length;

            if overlap > orf_len * 0.65
                && other_cod > gcod + 0.5
                && other_len as f64 > glen as f64 * 1.2
            {
                genes[i].shadow_pen = genes[i].shadow_pen.min(0.15);
                break;
            } else if overlap > orf_len * 0.45
                && other_cod > gcod + 0.8
                && other_len as f64 > glen as f64 * 1.8
            {
                genes[i].shadow_pen = genes[i].shadow_pen.min(0.35);
            }
        }
    }
}

pub fn operon_boost(genes: &mut [Gene], window: usize, min_score: f64) {
    let n = genes.len();
    if n < 2 { return; }
    let mut sorted: Vec<usize> = (0..n).collect();
    sorted.sort_by_key(|&i| genes[i].start);

    let mut boosts = vec![0.0f64; n];

    for ii in 0..sorted.len() {
        let i = sorted[ii];
        let mut nb = 0.0;

        for jj in (ii+1)..sorted.len() {
            let j = sorted[jj];
            if genes[j].start > genes[i].end + window { break; }
            if genes[j].is_plus == genes[i].is_plus && genes[j].score > min_score {
                let gap = if genes[j].start > genes[i].end { genes[j].start - genes[i].end } else { 0 };
                let df = (1.0 - gap as f64 / window as f64).max(0.0);
                nb += df * genes[j].score;
            }
        }

        for jj in (0..ii).rev() {
            let j = sorted[jj];
            if genes[i].start > genes[j].end + window { break; }
            if genes[j].is_plus == genes[i].is_plus && genes[j].score > min_score {
                let gap = if genes[i].start > genes[j].end { genes[i].start - genes[j].end } else { 0 };
                let df = (1.0 - gap as f64 / window as f64).max(0.0);
                nb += df * genes[j].score;
            }
        }

        if nb > 1.5 { boosts[i] = 0.08; }
        else if nb > 0.8 { boosts[i] = 0.05; }
        else if nb > 0.4 { boosts[i] = 0.03; }
    }

    for i in 0..n {
        genes[i].score += boosts[i];
    }
}

/// Intergenic distance model (based on Prodigal's intergenic_mod).
/// Models gene spacing explicitly in the DP:
/// - Same-strand, close (<60bp): bonus +0.04 to +0.08
/// - Same-strand, medium (60-180bp): neutral
/// - Same-strand, far (>180bp) or opposite strand: penalty -0.04
/// - Adjacent genes (0-4bp): cancel any RBS penalty (operon-internal)
/// Applied to DP weight so spacing affects gene selection.
pub fn connection_score(genes: &mut [Gene]) {
    let n = genes.len();
    if n < 2 { return; }

    let mut sorted: Vec<usize> = (0..n).collect();
    sorted.sort_by_key(|&i| genes[i].start);

    for ii in 0..sorted.len() {
        let i = sorted[ii];
        let mut best_conn = 0.0f64;

        // Check downstream neighbor (same strand)
        for jj in (ii + 1)..sorted.len() {
            let j = sorted[jj];
            if genes[j].start > genes[i].end + 500 { break; }
            if genes[j].is_plus != genes[i].is_plus { continue; }
            if genes[j].score < 0.30 { continue; }

            let gap = if genes[j].start > genes[i].end {
                genes[j].start - genes[i].end
            } else {
                0
            };
            let overlap = if genes[i].end >= genes[j].start {
                genes[i].end - genes[j].start + 1
            } else {
                0
            };

            let conn = if overlap >= 1 && overlap <= 4 {
                0.04  // ATGA-like coupling
            } else if overlap > 4 && overlap <= 30 {
                0.02
            } else if gap <= 20 {
                0.03  // very close
            } else if gap <= 50 {
                0.02
            } else if gap <= 150 {
                0.01
            } else {
                0.0
            };

            if conn > best_conn { best_conn = conn; }
            break;
        }

        // Check upstream neighbor (same strand)
        for jj in (0..ii).rev() {
            let j = sorted[jj];
            if genes[i].start > genes[j].end + 500 { break; }
            if genes[j].is_plus != genes[i].is_plus { continue; }
            if genes[j].score < 0.30 { continue; }

            let gap = if genes[i].start > genes[j].end {
                genes[i].start - genes[j].end
            } else {
                0
            };
            let overlap = if genes[j].end >= genes[i].start {
                genes[j].end - genes[i].start + 1
            } else {
                0
            };

            let conn = if overlap >= 1 && overlap <= 4 {
                0.04
            } else if overlap > 4 && overlap <= 30 {
                0.02
            } else if gap <= 20 {
                0.03
            } else if gap <= 50 {
                0.02
            } else if gap <= 150 {
                0.01
            } else {
                0.0
            };

            if conn > best_conn { best_conn = conn; }
            break;
        }

        genes[i].adj_bonus = best_conn;
        genes[i].weight += best_conn;
    }
}

/// Operon-aware rescue: if a gap exists between two predicted same-strand genes,
/// and an ORF fills that gap, it's likely part of the operon.
/// Operons share one promoter — internal genes may lack their own RBS.
/// This rescues genes that fail RBS threshold but sit inside operons.
pub fn operon_rescue(selected: &mut Vec<usize>, all: &[Gene]) {
    if selected.len() < 3 { return; }

    // Sort selected by position
    let mut sel_sorted: Vec<(usize, usize, usize, bool)> = selected.iter()
        .map(|&i| (all[i].start, all[i].end, i, all[i].is_plus))
        .collect();
    sel_sorted.sort();

    let sel_set: std::collections::HashSet<usize> = selected.iter().copied().collect();

    let mut rescued = Vec::new();

    // Find gaps between same-strand consecutive genes
    for w in sel_sorted.windows(2) {
        let (_, end1, _, strand1) = w[0];
        let (start2, _, _, strand2) = w[1];
        if strand1 != strand2 { continue; }

        let gap = if start2 > end1 { start2 - end1 } else { continue };
        if gap < 100 || gap > 3000 { continue; } // reasonable operon gap

        // Look for ORFs that fill this gap
        for (i, orf) in all.iter().enumerate() {
            if sel_set.contains(&i) { continue; }
            if orf.is_plus != strand1 { continue; }
            if orf.start < end1.saturating_sub(4) { continue; }
            if orf.end > start2 + 4 { continue; }
            if orf.length < 90 { continue; }

            // Moderate criteria: operon internal genes still need decent coding signal
            // But RBS requirement is relaxed (no own promoter needed)
            let min_score = 0.35;
            if orf.score < min_score { continue; }
            if !orf.is_longest { continue; }
            if !orf.is_atg() && orf.length < 200 { continue; } // non-ATG short = risky

            // Check no excessive overlap with neighbors
            let ov1 = if orf.start < end1 { end1 - orf.start + 1 } else { 0 };
            let ov2 = if orf.end > start2 { orf.end - start2 + 1 } else { 0 };
            if ov1 > 50 || ov2 > 50 { continue; }

            rescued.push((i, orf.score));
        }
    }

    // Sort by score, add top candidates
    rescued.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    rescued.dedup_by_key(|x| x.0);
    for &(idx, _) in rescued.iter().take(100) {
        selected.push(idx);
    }
}

// keep_best_starts: REMOVED (dead code, never called).
// Active version is inline in pipeline.rs run_prediction() step 2.

pub fn rescue_atypical(selected: &mut Vec<usize>, all: &[Gene]) {
    let sel_set: HashSet<(usize, usize, bool)> = selected.iter().map(|&i| {
        (all[i].start, all[i].end, all[i].is_plus)
    }).collect();
    let mut candidates: Vec<(usize, f64)> = Vec::new();

    for (i, orf) in all.iter().enumerate() {
        if sel_set.contains(&(orf.start, orf.end, orf.is_plus)) { continue; }
        if orf.shadow_pen < 0.5 { continue; }

        let rbs_pwm_norm = ((orf.rbs_pwm + 3.0) / 13.0).clamp(0.0, 1.0);
        let rbs_c = orf.rbs.max(rbs_pwm_norm * 0.85);

        let mut rescue = false;
        if rbs_c >= 0.55 && orf.length >= 180 && orf.is_longest && orf.is_atg() {
            rescue = true;
        } else if rbs_c >= 0.70 && orf.length >= 120 && orf.is_longest {
            rescue = true;
        } else if orf.hex_avg > 0.8 && orf.length >= 200 && orf.is_longest && orf.frame_bias > 0.5 {
            rescue = true;
        }

        if !rescue { continue; }

        let glen = (orf.end - orf.start + 1) as f64;
        let mut fits = true;
        for &si in selected.iter() {
            let a = &all[si];
            let ov_s = orf.start.max(a.start);
            let ov_e = orf.end.min(a.end);
            if ov_e >= ov_s {
                let overlap = (ov_e - ov_s + 1) as f64;
                if overlap > 30.0 || overlap / glen > 0.15 {
                    fits = false;
                    break;
                }
            }
        }

        if fits {
            let priority = orf.rbs * 0.4 + (orf.length as f64 / 3000.0) * 0.3 + orf.score * 0.3;
            candidates.push((i, priority));
        }
    }

    candidates.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    for &(i, _) in candidates.iter().take(60) {
        selected.push(i);
    }
}

// refine_starts + start_quality: REMOVED (dead code, disabled since F1 0.927→0.833).
// Start ranking is now done by LightGBM V2 model (pipeline.rs step after DP).
