use std::collections::HashMap;
use std::fs;

use crate::types::Gene;
use crate::io::rev_comp;
use crate::orf::find_orfs;
use crate::coding::{
    train_hex_initial, train_hex_from_set, score_hex_all, blend_hex, merge_hex,
    train_intergenic_hex, train_mono, train_mono_from_set, score_mono_all,
    merge_mono, blend_mono, train_atypical_hex, score_atypical,
};
use crate::rbs::{score_rbs_at, build_rbs_pwm, score_rbs_pwm, compute_upstream_at, score_leaderless, score_start_context};
use crate::scoring::{composite_score, compute_gc3, compute_gc3_bias, edge_coding_score};
use crate::selection::{
    dp_select, gap_fill, detect_shadows, connection_score, operon_rescue,
    rescue_atypical, refine_starts, gap_targeted_rescue, filter_same_strand_overlaps,
};
use crate::types::{HexModel, MonoModel};

pub fn main_cli() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: prokrustes <fasta> [--ncrna <gff>] [--debug <start> <end> <strand>] [--dump-all-orfs] [--dump-tsv] [--hmm]");
        std::process::exit(1);
    }
    let genome = crate::io::read_fasta(&args[1]);

    // Parse --ncrna flag: GFF file with rRNA/tRNA regions to mask
    let ncrna_regions = {
        let mut regions = Vec::new();
        let mut i = 2;
        while i < args.len() {
            if args[i] == "--ncrna" && i + 1 < args.len() {
                regions = parse_ncrna_gff(&args[i + 1]);
                break;
            }
            i += 1;
        }
        regions
    };

    // Terminator detection mode
    if args.iter().any(|a| a == "--terminators") {
        let rc_seq = crate::io::rev_comp(&genome);
        let terms = crate::terminator::find_all_terminators(&genome, &rc_seq);
        eprintln!("Found {} Rho-independent terminators", terms.len());
        println!("##gff-version 3");
        for t in &terms {
            let strand = if t.strand { "+" } else { "-" };
            println!(".\tterminator_finder\tterminator\t{}\t{}\t{:.1}\t{}\t.\thairpin={:.1};tail={:.1};conf={:.0}",
                t.start, t.end, -t.hairpin_score, strand, t.hairpin_score, t.tail_score, t.confidence);
        }
        return;
    }

    // Debug mode: trace a specific gene region
    if args.len() >= 5 && args[2] == "--debug" {
        let target_start: usize = args[3].parse().unwrap_or(0);
        let target_end: usize = args[4].parse().unwrap_or(0);
        let target_strand = if args.len() > 5 { &args[5] } else { "+" };
        debug_gene(&genome, target_start, target_end, target_strand == "+");
        return;
    }

    // HMM mode: full Viterbi gene finder
    if args.iter().any(|a| a == "--hmm") {
        let rc_seq = crate::io::rev_comp(&genome);
        let plus_orfs = crate::orf::find_orfs(&genome, true, genome.len());
        let minus_orfs = crate::orf::find_orfs(&rc_seq, false, genome.len());
        let hex = {
            let t1p = crate::coding::train_hex_initial(&genome, &plus_orfs, 900, true);
            let t1m = crate::coding::train_hex_initial(&rc_seq, &minus_orfs, 900, true);
            crate::coding::merge_hex(&t1p, &t1m).unwrap_or([0.0; crate::types::N_HEX])
        };
        let mut genes = crate::viterbi::hmm_gene_finder(&genome, &rc_seq, &hex, &hex, None);
        filter_ncrna_overlaps(&mut genes, &ncrna_regions);
        print!("[");
        for (i, g) in genes.iter().enumerate() {
            if i > 0 { print!(","); }
            print!("{{\"start\":{},\"end\":{},\"strand\":\"{}\"}}", g.start, g.end, if g.is_plus { "+" } else { "-" });
        }
        println!("]");
        return;
    }

    // Dump TSV: output ALL scored ORFs with all features
    if args.len() >= 3 && args[2] == "--dump-all-orfs" {
        let all = annotate_return_all(&genome);
        dump_tsv(&all);
        return;
    }

    // Dump TSV: output only final selected genes with all features
    if args.len() >= 3 && args[2] == "--dump-tsv" {
        let (genes, _) = annotate(&genome);
        dump_tsv(&genes);
        return;
    }

    let (mut genes, all_orfs) = annotate(&genome);
    filter_ncrna_overlaps(&mut genes, &ncrna_regions);

    // IS element / repeat family correction (post-annotation)
    crate::is_elements::correct_repeat_families(&genome, &mut genes);

    // Parse --model flag or auto-detect LightGBM start ranker
    let model_path = {
        let mut path = None;
        let mut i = 2;
        while i < args.len() {
            if args[i] == "--model" && i + 1 < args.len() {
                path = Some(args[i + 1].clone());
                break;
            }
            i += 1;
        }
        // Auto-detect model if not specified
        if path.is_none() {
            let candidates = [
                "models/start_ranker.lgb",
                "start_ranker.lgb",
                "/opt/prokrustes/models/start_ranker.lgb",
                "../models/start_ranker.lgb",
            ];
            // Also try relative to the binary location
            let mut all_candidates: Vec<String> = candidates.iter().map(|s| s.to_string()).collect();
            if let Ok(exe) = std::env::current_exe() {
                if let Some(dir) = exe.parent() {
                    all_candidates.push(format!("{}/models/start_ranker.lgb", dir.display()));
                    all_candidates.push(format!("{}/../models/start_ranker.lgb", dir.display()));
                }
            }
            for candidate in &all_candidates {
                if std::path::Path::new(candidate).exists() {
                    eprintln!("Auto-detected model: {}", candidate);
                    path = Some(candidate.clone());
                    break;
                }
            }
        }
        path
    };

    // LGB start re-ranking: use scored ORFs from dump-all-orfs
    // The model needs fully computed features (rbs, hex_avg, score, etc.)
    // So we get them via annotate_return_all() which scores everything
    if let Some(ref model_path) = model_path {
        if let Some(model) = crate::lgb_model::LgbModel::load(model_path) {
            // Index all scored ORFs by (strand, end) = same stop group
            let mut end_idx: std::collections::HashMap<(bool, usize), Vec<&crate::types::Gene>> =
                std::collections::HashMap::new();
            for orf in &all_orfs {
                end_idx.entry((orf.is_plus, orf.end)).or_default().push(orf);
            }

            let mut corrected = 0;
            for gene in genes.iter_mut() {
                let candidates = match end_idx.get(&(gene.is_plus, gene.end)) {
                    Some(c) if c.len() >= 2 => c,
                    _ => continue,
                };

                let max_len = candidates.iter().map(|c| c.length).max().unwrap_or(1);
                let max_rbs = candidates.iter().map(|c| c.rbs).fold(0.0f64, f64::max);
                let max_hex = candidates.iter().map(|c| c.hex_avg).fold(f64::NEG_INFINITY, f64::max);
                let max_score = candidates.iter().map(|c| c.score).fold(0.0f64, f64::max);
                let n_cands = candidates.len() as f64;

                let mut best_lgb = f64::NEG_INFINITY;
                let mut best_c: Option<&&crate::types::Gene> = None;

                for (ci, c) in candidates.iter().enumerate() {
                    let features = [
                        c.rbs, c.rbs_pwm, c.start_type(), c.hex_avg, c.frame_bias, c.edge,
                        c.length as f64, c.score, c.start_nn, c.start_ctx,
                        c.upstream_at, c.leaderless, c.gc3_bias,
                        if c.is_longest { 1.0 } else { 0.0 },
                        ci as f64, n_cands,
                        c.rbs - max_rbs,
                        c.hex_avg - max_hex,
                        (c.length as f64) - (max_len as f64),
                        c.score - max_score,
                        c.length as f64 / max_len as f64,
                        c.upstream_coding,
                        (max_len - c.length) as f64 / 3.0,
                        1.0 - c.length as f64 / max_len as f64,
                    ];
                    let s = model.predict(&features);
                    if s > best_lgb {
                        best_lgb = s;
                        best_c = Some(c);
                    }
                }

                if let Some(best) = best_c {
                    if best.start != gene.start {
                        gene.start = best.start;
                        gene.end = best.end;
                        gene.length = best.length;
                        gene.seq_start = best.seq_start;
                        gene.seq_end = best.seq_end;
                        gene.start_codon = best.start_codon;
                        corrected += 1;
                    }
                }
            }
            if corrected > 0 {
                eprintln!("LGB start ranker: {} starts corrected", corrected);
            }
        }
    }

    // Parse --compare flag: comparison genomes for conservation scoring
    let compare_fastas: Vec<String> = {
        let mut fastas = Vec::new();
        let mut i = 2;
        while i < args.len() {
            if args[i] == "--compare" {
                i += 1;
                while i < args.len() && !args[i].starts_with("--") {
                    fastas.push(args[i].clone());
                    i += 1;
                }
                break;
            }
            i += 1;
        }
        fastas
    };

    // Conservation scoring + start correction
    let cons_scores = if !compare_fastas.is_empty() {
        eprintln!("Comparative genomics: {} comparison genomes", compare_fastas.len());
        let mut comp_data: Vec<(Vec<u8>, Vec<Gene>, crate::conservation::KmerIndex)> = Vec::new();
        for fasta in &compare_fastas {
            let comp_genome = crate::io::read_fasta(fasta);
            let (comp_genes, _) = annotate(&comp_genome);
            let idx = crate::conservation::KmerIndex::build(&comp_genome, &comp_genes);
            eprintln!("  {}: {} genes", fasta.rsplit('/').next().unwrap_or(fasta), comp_genes.len());
            comp_data.push((comp_genome, comp_genes, idx));
        }
        let idx_refs: Vec<&crate::conservation::KmerIndex> = comp_data.iter().map(|(_, _, idx)| idx).collect();

        // Score conservation
        let scores = crate::conservation::score_conservation(&genome, &genes, &idx_refs);
        let conserved = scores.iter().filter(|&&s| s > 0.3).count();
        let unique = scores.iter().filter(|&&s| s <= 0.1).count();
        eprintln!("  Conserved: {}, Unique: {}", conserved, unique);

        // Start correction: extend truncated genes to match ortholog length
        let comp_refs: Vec<(&[u8], &[Gene])> = comp_data.iter()
            .map(|(g, genes, _)| (g.as_slice(), genes.as_slice())).collect();
        // Need all_orfs for finding alternative starts in same stop group
        let mut all_orfs_for_correction: Vec<Gene> = Vec::new();
        all_orfs_for_correction.extend(find_orfs(&genome, true, genome.len()).into_iter());
        all_orfs_for_correction.extend(find_orfs(&crate::io::rev_comp(&genome), false, genome.len()).into_iter());

        let n_corrected = crate::conservation::correct_starts_by_orthologs(
            &genome, &mut genes, &all_orfs_for_correction, &idx_refs, &comp_refs);
        if n_corrected > 0 {
            eprintln!("  Start correction by orthologs: {} genes extended", n_corrected);
        }

        Some(scores)
    } else {
        None
    };

    // Output GFF3
    println!("##gff-version 3");
    for (i, g) in genes.iter().enumerate() {
        let strand = if g.is_plus { "+" } else { "-" };
        let codon = std::str::from_utf8(&g.start_codon).unwrap_or("???");
        let mut attrs = format!(
            "ID=gene_{};start_type={};rbs_score={:.2};hex_score={:.4}",
            i + 1, codon, g.rbs, g.hex_avg
        );
        if let Some(ref scores) = cons_scores {
            if i < scores.len() {
                let label = crate::conservation::conservation_label(scores[i]);
                attrs.push_str(&format!(";conservation={}", label));
            }
        }
        println!(".\tprokrustes\tCDS\t{}\t{}\t{:.1}\t{}\t0\t{}",
            g.start, g.end, g.score, strand, attrs);
    }
}

/// Dump all genes as TSV with every internal feature.
fn dump_tsv(genes: &[Gene]) {
    println!("start\tend\tstrand\tlength\tstart_codon\tis_longest\thex_avg\thex_total\tframe_bias\thex_cov\trbs\trbs_pwm\tgc3\tupstream_at\tmono\tedge\tscore\tweight\tshadow_pen\tleaderless\tstart_ctx\tgc3_bias\tadj_bonus\tviterbi_frac\tstart_nn\tstop_group\tframe");
    for g in genes {
        let codon = std::str::from_utf8(&g.start_codon).unwrap_or("???");
        let strand = if g.is_plus { "+" } else { "-" };
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{}\t{}",
            g.start, g.end, strand, g.length, codon, g.is_longest,
            g.hex_avg, g.hex_total, g.frame_bias, g.hex_cov,
            g.rbs, g.rbs_pwm, g.gc3, g.upstream_at,
            g.mono, g.edge, g.score, g.weight, g.shadow_pen,
            g.leaderless, g.start_ctx, g.gc3_bias, g.adj_bonus, g.viterbi_frac, g.start_nn,
            g.stop_group, g.frame);
    }
}

/// Run full annotation pipeline but return ALL scored ORFs (not just selected).
pub fn annotate_return_all(genome: &[u8]) -> Vec<Gene> {
    let glen = genome.len();
    let rc = crate::io::rev_comp(genome);
    let mut plus = crate::orf::find_orfs(genome, true, glen);
    let mut minus = crate::orf::find_orfs(&rc, false, glen);

    // Same scoring pipeline as annotate() but skip selection
    for orf in plus.iter_mut() {
        orf.gc3 = crate::scoring::compute_gc3(genome, orf.seq_start, orf.seq_end);
        orf.gc3_bias = crate::scoring::compute_gc3_bias(genome, orf.seq_start, orf.seq_end);
        orf.rbs = crate::rbs::score_rbs_at(genome, orf.seq_start);
        orf.upstream_at = crate::rbs::compute_upstream_at(genome, orf.seq_start);
        orf.leaderless = crate::rbs::score_leaderless(genome, orf.seq_start);
        orf.start_ctx = crate::rbs::score_start_context(genome, orf.seq_start);
    }
    for orf in minus.iter_mut() {
        orf.gc3 = crate::scoring::compute_gc3(&rc, orf.seq_start, orf.seq_end);
        orf.gc3_bias = crate::scoring::compute_gc3_bias(&rc, orf.seq_start, orf.seq_end);
        orf.rbs = crate::rbs::score_rbs_at(&rc, orf.seq_start);
        orf.upstream_at = crate::rbs::compute_upstream_at(&rc, orf.seq_start);
        orf.leaderless = crate::rbs::score_leaderless(&rc, orf.seq_start);
        orf.start_ctx = crate::rbs::score_start_context(&rc, orf.seq_start);
    }

    // Quick hexamer training on long ORFs
    let t1p = train_hex_initial(genome, &plus, 900, true);
    let t1m = train_hex_initial(&rc, &minus, 900, true);
    if let Some(hex) = merge_hex(&t1p, &t1m) {
        let gc3_target = 0.5;
        score_hex_all(genome, &mut plus, &hex);
        score_hex_all(&rc, &mut minus, &hex);
        for orf in plus.iter_mut() {
            orf.score = composite_score(orf, gc3_target);
        }
        for orf in minus.iter_mut() {
            orf.score = composite_score(orf, gc3_target);
        }
    }

    // Merge and return all
    let mut all: Vec<Gene> = Vec::new();
    all.extend(plus);
    all.extend(minus);
    all.sort_by_key(|g| g.start);
    all
}

/// Debug mode: trace WHY a specific gene is found or missed.
fn debug_gene(genome: &[u8], target_start: usize, target_end: usize, target_plus: bool) {
    let glen = genome.len();
    let rc = crate::io::rev_comp(genome);
    let strand_seq = if target_plus { genome } else { &rc };
    let strand_str = if target_plus { "+" } else { "-" };

    eprintln!("═══ DEBUG: gene {}..{} ({}) ═══", target_start, target_end, strand_str);
    eprintln!("  Length: {}bp = {} aa", target_end - target_start + 1, (target_end - target_start - 2) / 3);
    eprintln!();

    // 1. Find ORFs in this region
    let orfs = crate::orf::find_orfs(strand_seq, target_plus, glen);
    let nearby: Vec<&Gene> = orfs.iter().filter(|o| {
        o.start <= target_end + 100 && o.end >= target_start.saturating_sub(100)
    }).collect();

    eprintln!("  ORFs near target region: {}", nearby.len());
    let exact = nearby.iter().filter(|o| {
        (o.start as i64 - target_start as i64).unsigned_abs() < 50 &&
        (o.end as i64 - target_end as i64).unsigned_abs() < 50
    }).collect::<Vec<_>>();

    if exact.is_empty() {
        eprintln!("  ✗ NO ORF matches target coordinates!");
        eprintln!("  Nearby ORFs:");
        for o in nearby.iter().take(10) {
            let seq_start_codon = if o.seq_start + 3 <= strand_seq.len() {
                std::str::from_utf8(&strand_seq[o.seq_start..o.seq_start+3]).unwrap_or("???")
            } else { "???" };
            eprintln!("    {}..{} ({}) {}bp start={} longest={}",
                o.start, o.end, strand_str, o.length, seq_start_codon, o.is_longest);
        }
    } else {
        eprintln!("  ✓ Found matching ORF(s):");
        for o in &exact {
            let seq_start_codon = if o.seq_start + 3 <= strand_seq.len() {
                std::str::from_utf8(&strand_seq[o.seq_start..o.seq_start+3]).unwrap_or("???")
            } else { "???" };
            eprintln!("    {}..{} {}bp start={} longest={}", o.start, o.end, o.length, seq_start_codon, o.is_longest);
        }
    }

    // 1b. Show scores of matching ORFs (after full pipeline training)
    eprintln!();
    eprintln!("  Pre-filter scores (from full pipeline):");
    {
        let rc = crate::io::rev_comp(genome);
        let all_orfs = if target_plus {
            crate::orf::find_orfs(genome, true, glen)
        } else {
            crate::orf::find_orfs(&rc, false, glen)
        };
        for o in &all_orfs {
            if (o.start as i64 - target_start as i64).unsigned_abs() < 50 &&
               (o.end as i64 - target_end as i64).unsigned_abs() < 50 {
                let seq = if target_plus { genome } else { &rc };
                let rbs_val = crate::rbs::score_rbs_at(seq, o.seq_start);
                eprintln!("    {}..{} {}bp rbs={:.2} is_longest={} is_atg={}",
                    o.start, o.end, o.length, rbs_val, o.is_longest, o.is_atg());
            }
        }
    }

    // 2. Run full annotation and check what happens to this region
    eprintln!();
    eprintln!("  Running full annotation...");
    let (genes, _) = annotate(genome);

    let matching: Vec<&Gene> = genes.iter().filter(|g| {
        g.start <= target_end && g.end >= target_start && g.is_plus == target_plus
    }).collect();

    if matching.is_empty() {
        eprintln!("  ✗ GENE NOT IN FINAL OUTPUT!");
        // Check what's nearby in output
        let near_output: Vec<&Gene> = genes.iter().filter(|g| {
            g.start <= target_end + 500 && g.end >= target_start.saturating_sub(500)
        }).collect();
        eprintln!("  Nearby predictions (±500bp):");
        for g in near_output.iter().take(5) {
            let overlap_start = g.start.max(target_start);
            let overlap_end = g.end.min(target_end);
            let overlap = if overlap_end > overlap_start { overlap_end - overlap_start } else { 0 };
            eprintln!("    {}..{} ({}) {}bp score={:.3} overlap={}bp",
                g.start, g.end, if g.is_plus { "+" } else { "-" },
                g.end - g.start + 1, g.score, overlap);
        }
    } else {
        eprintln!("  ✓ GENE FOUND in output:");
        for g in &matching {
            let overlap_start = g.start.max(target_start);
            let overlap_end = g.end.min(target_end);
            let overlap = if overlap_end > overlap_start { overlap_end - overlap_start } else { 0 };
            let target_len = target_end - target_start + 1;
            let pct = overlap as f64 / target_len as f64 * 100.0;
            eprintln!("    {}..{} {}bp score={:.3} hex={:.3} rbs={:.2} fb={:.2} overlap={}bp ({:.0}%)",
                g.start, g.end, g.end - g.start + 1, g.score, g.hex_avg, g.rbs, g.frame_bias, overlap, pct);
        }
    }
}

/// Score all ORFs: hexamer, monocodon, edge (optional), and composite.
fn score_all_orfs(
    plus: &mut [Gene], minus: &mut [Gene],
    genome: &[u8], rc: &[u8],
    hex_model: &HexModel, mono_model: &Option<MonoModel>,
    gc3_target: f64, compute_edge: bool,
) {
    score_hex_all(genome, plus, hex_model);
    score_hex_all(rc, minus, hex_model);

    if let Some(mt) = mono_model {
        score_mono_all(genome, plus, mt);
        score_mono_all(rc, minus, mt);
    }

    if compute_edge {
        for orf in plus.iter_mut() {
            orf.edge = edge_coding_score(genome, orf.seq_start, hex_model, 36);
        }
        for orf in minus.iter_mut() {
            orf.edge = edge_coding_score(rc, orf.seq_start, hex_model, 36);
        }
    }

    for orf in plus.iter_mut() {
        orf.score = composite_score(orf, gc3_target);
    }
    for orf in minus.iter_mut() {
        orf.score = composite_score(orf, gc3_target);
    }
}

/// Run prediction pipeline: filter → keep_best_starts → operon_boost → weights → DP → gap_fill.
/// Matches monolith's run_prediction exactly.
fn run_prediction(all_orfs: &mut Vec<Gene>, thresh_adj: f64) -> (Vec<usize>, Vec<usize>) {
    // 1. Filter by various criteria
    let mut filtered: Vec<usize> = Vec::new();
    for (i, orf) in all_orfs.iter().enumerate() {
        let length = orf.length;
        let cod = orf.hex_avg;
        let rbs = orf.rbs;
        let fb = orf.frame_bias;

        if cod < -0.3 && length < 600 { continue; }
        if length < 250 && cod < 0.05 && rbs < 0.15 { continue; }
        if fb < -0.4 && length < 500 { continue; }
        // Short gene filter: allow if RBS is strong (rpmJ=0.60 was killed here)
        if length < 180 && fb < 0.10 && cod < 0.20 && rbs < 0.40 { continue; }

        let mut thresh = if length >= 900 { 0.10 }
            else if length >= 600 { 0.19 }
            else if length >= 300 { 0.31 }
            else if length >= 150 { 0.39 }
            else { 0.49 };
        thresh += thresh_adj;

        if rbs > 0.60 { thresh *= 0.80; }
        else if rbs > 0.35 { thresh *= 0.90; }

        if orf.score >= thresh {
            filtered.push(i);
        }
    }

    // 2. Keep best 4 starts per stop group
    {
        let mut groups: HashMap<(bool, u32), Vec<usize>> = HashMap::new();
        for &i in &filtered {
            let orf = &all_orfs[i];
            let key = (orf.is_plus, orf.stop_group);
            groups.entry(key).or_default().push(i);
        }
        let mut selected_indices = Vec::new();
        for (_, indices) in &groups {
            let mut scored: Vec<(usize, f64)> = indices.iter().map(|&i| {
                let orf = &all_orfs[i];
                let rbs_pwm_norm = ((orf.rbs_pwm + 3.0) / 13.0).clamp(0.0, 1.0);
                let rbs_c = orf.rbs.max(rbs_pwm_norm * 0.85);
                let sc = orf.start_type();
                let cn = ((orf.hex_avg + 1.0) / 4.0).clamp(0.0, 1.0);
                let fb = (orf.frame_bias / 3.0).max(0.0);
                let edge = ((orf.edge + 0.5) / 3.0).clamp(0.0, 1.0);
                // Truncation penalty: log(lost_codons).
                // Frame-aware version (upstream_coding as feature) attempted but
                // worsened F1 (-7 TP). Need pairwise ranking, not additive penalty.
                // See TASK-006: reformulate as group-wise ranking within stop group.
                let max_len_in_group = indices.iter()
                    .map(|&j| all_orfs[j].length)
                    .max().unwrap_or(orf.length);
                let truncation_penalty = if orf.length < max_len_in_group {
                    let lost_codons = (max_len_in_group - orf.length) / 3;
                    0.015 * (lost_codons as f64 + 1.0).ln()
                } else {
                    0.0
                };
                let longest_bonus = if orf.is_longest { 0.06 } else { 0.0 };
                let nn_bonus = if orf.start_nn > 0.5 { (orf.start_nn - 0.5) * 0.40 } else { -0.03 };
                // N-terminal protein composition (stored in stop_nn field)
                // Positive = looks like real N-terminus, negative = looks like mid-protein
                let nterm_bonus = (orf.stop_nn * 0.08).clamp(-0.04, 0.08);
                let sq = 0.18 * rbs_c + 0.09 * cn + 0.06 * fb + 0.06 * edge
                    + 0.14 * sc + 0.12 * (1.0 - (-(orf.length as f64) / 500.0).exp())
                    + nn_bonus + nterm_bonus + longest_bonus - truncation_penalty;
                (i, sq)
            }).collect();
            scored.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
            for &(idx, sq) in scored.iter().take(4) {
                all_orfs[idx]._sq = sq;
                selected_indices.push(idx);
            }
        }
        filtered = selected_indices;
    }

    // 3. Operon boost on filtered subset
    {
        let mut sorted_idx = filtered.clone();
        sorted_idx.sort_by_key(|&i| all_orfs[i].start);
        let n = sorted_idx.len();
        let mut boosts = vec![0.0f64; n];
        for ii in 0..n {
            let i = sorted_idx[ii];
            let mut nb_score = 0.0;
            for jj in (ii + 1)..n {
                let j = sorted_idx[jj];
                if all_orfs[j].start > all_orfs[i].end + 250 { break; }
                let gap = if all_orfs[j].start > all_orfs[i].end {
                    all_orfs[j].start - all_orfs[i].end
                } else { 0 };
                if all_orfs[j].is_plus == all_orfs[i].is_plus && all_orfs[j].score > 0.42 {
                    let df = (1.0 - gap as f64 / 250.0).max(0.0);
                    nb_score += df * all_orfs[j].score;
                }
            }
            for jj in (0..ii).rev() {
                let j = sorted_idx[jj];
                if all_orfs[i].start > all_orfs[j].end + 250 { break; }
                let gap = if all_orfs[i].start > all_orfs[j].end {
                    all_orfs[i].start - all_orfs[j].end
                } else { 0 };
                if all_orfs[j].is_plus == all_orfs[i].is_plus && all_orfs[j].score > 0.42 {
                    let df = (1.0 - gap as f64 / 250.0).max(0.0);
                    nb_score += df * all_orfs[j].score;
                }
            }
            if nb_score > 1.5 { boosts[ii] = 0.08; }
            else if nb_score > 0.8 { boosts[ii] = 0.05; }
            else if nb_score > 0.4 { boosts[ii] = 0.03; }
        }
        for (ii, &idx) in sorted_idx.iter().enumerate() {
            all_orfs[idx].score += boosts[ii];
        }
    }

    // 3b. Connection scoring — DISABLED: F1 0.938→0.937
    // Doesn't improve over ML-weighted baseline. Borderline genes that
    // connection scoring helps are already captured by operon_boost.
    // Code kept in selection.rs for future use.
    // connection_score(&mut filt_genes);

    // 4. Compute weights (matches monolith exactly)
    for &i in &filtered {
        let orf = &all_orfs[i];
        let length = orf.length;
        let mut base = if length >= 900 { 0.12 }
            else if length >= 600 { 0.21 }
            else if length >= 300 { 0.33 }
            else if length >= 150 { 0.41 }
            else { 0.51 };
        base += thresh_adj;
        if orf.rbs > 0.60 { base *= 0.80; }
        else if orf.rbs > 0.35 { base *= 0.90; }

        let w = (orf.score - base) + (orf.hex_total * 0.004).max(0.0);
        all_orfs[i].weight = w;
    }

    // 5. Per-strand iterative DP (4 iterations)
    let plus_sel: Vec<usize> = filtered.iter().filter(|&&i| all_orfs[i].is_plus).copied().collect();
    let minus_sel: Vec<usize> = filtered.iter().filter(|&&i| !all_orfs[i].is_plus).copied().collect();

    let run_iter_dp = |indices: &[usize], all: &mut Vec<Gene>| -> Vec<usize> {
        if indices.is_empty() { return vec![]; }
        let mut sub: Vec<Gene> = indices.iter().map(|&i| all[i].clone()).collect();
        for orf in sub.iter_mut() {
            orf._base_weight = orf.weight;
            orf.adj_bonus = 0.0;
        }
        let mut results = Vec::new();
        for it in 0..4 {
            for orf in sub.iter_mut() {
                orf.weight = orf._base_weight + orf.adj_bonus;
            }
            results = dp_select(&sub, 60);
            if it < 3 {
                let mut sel_sorted: Vec<usize> = results.clone();
                sel_sorted.sort_by_key(|&i| sub[i].start);
                let sel_starts: Vec<usize> = sel_sorted.iter().map(|&i| sub[i].start).collect();
                for ci in 0..sub.len() {
                    let mut bonus = 0.0f64;
                    let target_lo = sub[ci].start.saturating_sub(300);
                    let idx = sel_starts.partition_point(|&s| s < target_lo);
                    for k in idx..sel_sorted.len().min(idx + 30) {
                        let si = sel_sorted[k];
                        if sub[si].start > sub[ci].end + 300 { break; }
                        if sub[si].is_plus != sub[ci].is_plus { continue; }
                        if sub[si].start == sub[ci].start && sub[si].end == sub[ci].end { continue; }
                        let gap: isize;
                        if sub[si].end <= sub[ci].start {
                            gap = (sub[ci].start - sub[si].end) as isize;
                        } else if sub[ci].end <= sub[si].start {
                            gap = (sub[si].start - sub[ci].end) as isize;
                        } else {
                            let ov = sub[si].end.min(sub[ci].end) as isize
                                - sub[si].start.max(sub[ci].start) as isize + 1;
                            if ov <= 4 { gap = -ov; } else { continue; }
                        }
                        let b = if gap >= -4 && gap <= 4 { 0.15 }
                            else if gap <= 30 { 0.10 }
                            else if gap <= 80 { 0.05 }
                            else if gap <= 150 { 0.02 }
                            else { 0.0 };
                        bonus = bonus.max(b);
                    }
                    sub[ci].adj_bonus = bonus;
                }
            }
        }
        // Copy back weights
        for (li, &oi) in indices.iter().enumerate() {
            all[oi].weight = sub[li].weight;
            all[oi].adj_bonus = sub[li].adj_bonus;
            all[oi]._base_weight = sub[li]._base_weight;
        }
        results.iter().map(|&si| indices[si]).collect()
    };

    let plus_res = run_iter_dp(&plus_sel, all_orfs);
    let minus_res = run_iter_dp(&minus_sel, all_orfs);

    // 6. Combined DP
    let combined: Vec<usize> = plus_res.iter().chain(minus_res.iter()).copied().collect();
    let combined_orfs: Vec<Gene> = combined.iter().map(|&i| all_orfs[i].clone()).collect();
    let final_sel = dp_select(&combined_orfs, 50);
    let mut results: Vec<usize> = final_sel.iter().map(|&si| combined[si]).collect();

    // 7. Gap fill
    gap_fill(&mut results, all_orfs, &filtered, 0.40, 120);

    (results, filtered)
}

pub fn annotate(genome: &[u8]) -> (Vec<Gene>, Vec<Gene>) {
    let glen = genome.len();
    let rc = rev_comp(genome);

    // 0. Auto-detect genetic code
    let code = crate::orf::detect_genetic_code(genome);
    match code {
        crate::orf::GeneticCode::Standard => eprintln!("Genetic code: 11 (standard bacterial)"),
        crate::orf::GeneticCode::Code4 => eprintln!("Genetic code: 4 (Mycoplasma — TGA=Trp, stops: TAA/TAG only)"),
    }

    // 1. Find all ORFs on both strands
    let mut plus = crate::orf::find_orfs_code(genome, true, glen, code);
    let mut minus = crate::orf::find_orfs_code(&rc, false, glen, code);

    // 2. Compute GC3 target from long ORFs (before hexamer training, matching monolith)
    for orf in plus.iter_mut() {
        orf.gc3 = compute_gc3(genome, orf.seq_start, orf.seq_end);
        orf.gc3_bias = compute_gc3_bias(genome, orf.seq_start, orf.seq_end);
    }
    for orf in minus.iter_mut() {
        orf.gc3 = compute_gc3(&rc, orf.seq_start, orf.seq_end);
        orf.gc3_bias = compute_gc3_bias(&rc, orf.seq_start, orf.seq_end);
    }

    let gc3_target = {
        let long: Vec<f64> = plus.iter().chain(minus.iter())
            .filter(|o| o.length >= 900 && o.is_atg() && o.is_longest)
            .map(|o| o.gc3).collect();
        if long.len() >= 50 {
            let limit = long.len().min(500);
            long[..limit].iter().sum::<f64>() / limit as f64
        } else {
            let alt: Vec<f64> = plus.iter().chain(minus.iter())
                .filter(|o| o.length >= 600 && o.is_longest)
                .map(|o| o.gc3).collect();
            let limit = alt.len().min(500);
            if limit > 0 { alt[..limit].iter().sum::<f64>() / limit as f64 } else { 0.5 }
        }
    };

    // 3. RBS and upstream features (before hexamer training, matching monolith)
    for orf in plus.iter_mut() {
        orf.rbs = score_rbs_at(genome, orf.seq_start);
        orf.upstream_at = compute_upstream_at(genome, orf.seq_start);
        orf.leaderless = score_leaderless(genome, orf.seq_start);
        orf.start_ctx = score_start_context(genome, orf.seq_start);
    }
    for orf in minus.iter_mut() {
        orf.rbs = score_rbs_at(&rc, orf.seq_start);
        orf.upstream_at = compute_upstream_at(&rc, orf.seq_start);
        orf.leaderless = score_leaderless(&rc, orf.seq_start);
        orf.start_ctx = score_start_context(&rc, orf.seq_start);
    }

    // 4. Initial hexamer tables (matching monolith's fallback chain)
    let t1p = train_hex_initial(genome, &plus, 900, true);
    let t1m = train_hex_initial(&rc, &minus, 900, true);
    let initial_hex = merge_hex(&t1p, &t1m).or_else(|| {
        let t2p = train_hex_initial(genome, &plus, 600, true);
        let t2m = train_hex_initial(&rc, &minus, 600, true);
        merge_hex(&t2p, &t2m)
    });

    let mut hex_model = match initial_hex {
        Some(m) => m,
        None => {
            // Emergency fallback
            let all: Vec<Gene> = plus.into_iter().chain(minus.into_iter()).collect();
            let mut result: Vec<Gene> = all.iter()
                .filter(|o| o.length >= 300 && o.is_longest && o.is_atg())
                .cloned().collect();
            result.sort_by_key(|g| g.start);
            return (result, all);
        }
    };

    // 5. Initial monocodon tables
    let mono_p = train_mono(genome, &plus, 900);
    let mono_m = train_mono(&rc, &minus, 900);
    let mut mono_model = merge_mono(&mono_p, &mono_m);

    // 6. Initial scoring (no edge)
    score_all_orfs(&mut plus, &mut minus, genome, &rc, &hex_model, &mono_model, gc3_target, false);

    // Filter obvious noise
    plus.retain(|o| o.score > 0.05);
    minus.retain(|o| o.score > 0.05);

    // 7. Iterative self-training (matches monolith: filter by score, 12 iterations)
    let mut prev_hash: Option<f64> = None;
    for iteration in 0..12 {
        let conf_thresh = 0.46 - iteration as f64 * 0.02;
        let min_train = (450usize).saturating_sub(iteration * 30).max(180);

        let conf_p: Vec<Gene> = plus.iter().filter(|o| {
            o.score > conf_thresh && o.length >= min_train && o.is_longest
        }).cloned().collect();
        let conf_m: Vec<Gene> = minus.iter().filter(|o| {
            o.score > conf_thresh && o.length >= min_train && o.is_longest
        }).cloned().collect();

        let t2p = train_hex_from_set(genome, &conf_p, 150);
        let t2m = train_hex_from_set(&rc, &conf_m, 150);

        if t2p.is_some() || t2m.is_some() {
            let t2 = merge_hex(&t2p, &t2m);
            if let Some(new_hex) = t2 {
                let blend_w = (0.35 + iteration as f64 * 0.05).min(0.70);
                hex_model = blend_hex(&hex_model, &new_hex, blend_w);

                // Update monocodon every 3rd iteration
                if iteration % 3 == 2 {
                    let mt_p = train_mono_from_set(genome, &conf_p, 200);
                    let mt_m = train_mono_from_set(&rc, &conf_m, 200);
                    let mt = merge_mono(&mt_p, &mt_m);
                    if let Some(new_mt) = mt {
                        mono_model = Some(match mono_model {
                            Some(existing) => blend_mono(&existing, &new_mt, 0.35),
                            None => new_mt,
                        });
                    }
                }

                score_all_orfs(&mut plus, &mut minus, genome, &rc, &hex_model, &mono_model, gc3_target, false);
            }
        }

        // Convergence check
        let sample: f64 = hex_model.iter().take(50).sum();
        if let Some(ph) = prev_hash {
            if (sample - ph).abs() < 0.1 { break; }
        }
        prev_hash = Some(sample);
    }

    // 8. First intergenic refinement
    {
        let conf_genes: Vec<Gene> = plus.iter().chain(minus.iter())
            .filter(|o| o.score > 0.48 && o.length >= 350 && o.is_longest)
            .cloned().collect();

        if conf_genes.len() > 200 {
            let intergenic = get_intergenic_regions(genome, &conf_genes, 50);
            if !intergenic.is_empty() {
                let conf_p: Vec<Gene> = plus.iter()
                    .filter(|o| o.score > 0.48 && o.length >= 300 && o.is_longest)
                    .cloned().collect();
                let conf_m: Vec<Gene> = minus.iter()
                    .filter(|o| o.score > 0.48 && o.length >= 300 && o.is_longest)
                    .cloned().collect();

                let it_p = train_intergenic_hex(genome, &conf_p, &intergenic);
                let rc_ig: Vec<Vec<u8>> = intergenic.iter().map(|s| rev_comp(s)).collect();
                let it_m = train_intergenic_hex(&rc, &conf_m, &rc_ig);
                let it_model = merge_hex(&it_p, &it_m);

                if let Some(igm) = it_model {
                    hex_model = blend_hex(&hex_model, &igm, 0.30);
                    score_all_orfs(&mut plus, &mut minus, genome, &rc, &hex_model, &mono_model, gc3_target, false);

                    // Extra refinement rounds
                    for extra in 0..3 {
                        let ct = 0.40 - extra as f64 * 0.02;
                        let cp: Vec<Gene> = plus.iter()
                            .filter(|o| o.score > ct && o.length >= 200 && o.is_longest)
                            .cloned().collect();
                        let cm: Vec<Gene> = minus.iter()
                            .filter(|o| o.score > ct && o.length >= 200 && o.is_longest)
                            .cloned().collect();
                        let t3p = train_hex_from_set(genome, &cp, 150);
                        let t3m = train_hex_from_set(&rc, &cm, 150);
                        if let Some(new_m) = merge_hex(&t3p, &t3m) {
                            hex_model = blend_hex(&hex_model, &new_m, 0.30 + extra as f64 * 0.05);
                            score_all_orfs(&mut plus, &mut minus, genome, &rc, &hex_model, &mono_model, gc3_target, false);
                        }
                    }
                }
            }
        }
    }

    // 9. Final scoring with edge
    score_all_orfs(&mut plus, &mut minus, genome, &rc, &hex_model, &mono_model, gc3_target, true);

    // 9x. Compute upstream in-frame coding potential (frame-aware truncation evidence)
    for orf in plus.iter_mut() {
        orf.upstream_coding = crate::scoring::upstream_inframe_coding(genome, orf.seq_start, &hex_model, 500);
    }
    for orf in minus.iter_mut() {
        orf.upstream_coding = crate::scoring::upstream_inframe_coding(&rc, orf.seq_start, &hex_model, 500);
    }

    // 9a. HMM as additional signal: run full Viterbi, compute overlap with each ORF
    {
        let hmm_genes = crate::viterbi::hmm_gene_finder(genome, &rc, &hex_model, &hex_model, None);
        // For each pipeline ORF, compute fraction overlapping with ANY HMM gene
        let set_hmm_frac = |orfs: &mut [Gene]| {
            for orf in orfs.iter_mut() {
                let mut best_frac = 0.0f64;
                for hg in &hmm_genes {
                    if hg.is_plus != orf.is_plus { continue; }
                    let ov_s = orf.start.max(hg.start);
                    let ov_e = orf.end.min(hg.end);
                    if ov_e > ov_s {
                        let frac = (ov_e - ov_s) as f64 / orf.length.max(1) as f64;
                        if frac > best_frac { best_frac = frac; }
                    }
                }
                orf.viterbi_frac = best_frac;
            }
        };
        set_hmm_frac(&mut plus);
        set_hmm_frac(&mut minus);
    }

    // 9b. Atypical (AT-rich) gene rescue: train separate hexamer on low-GC ORFs
    {
        let genome_gc = {
            let gc = genome.iter().filter(|&&c| c == b'G' || c == b'C').count();
            gc as f64 / genome.len() as f64
        };
        let atyp_p = train_atypical_hex(genome, &plus, genome_gc);
        let atyp_m = train_atypical_hex(&rc, &minus, genome_gc);
        if let Some(ref model) = atyp_p {
            score_atypical(genome, &mut plus, model, genome_gc);
        }
        if let Some(ref model) = atyp_m {
            score_atypical(&rc, &mut minus, model, genome_gc);
        }
        // Re-compute composite scores for boosted ORFs
        for orf in plus.iter_mut() {
            orf.score = composite_score(orf, gc3_target);
        }
        for orf in minus.iter_mut() {
            orf.score = composite_score(orf, gc3_target);
        }
    }

    // 9c. Start + stop context neural models
    {
        let sm_p = crate::start_model::train_start_model(genome, &plus);
        let sm_m = crate::start_model::train_start_model(&rc, &minus);
        if let Some(ref model) = sm_p {
            crate::start_model::score_start_model(genome, &mut plus, model);
        }
        if let Some(ref model) = sm_m {
            crate::start_model::score_start_model(&rc, &mut minus, model);
        }
        // N-terminal protein composition model
        let nt_p = crate::start_model::train_nterm_model(genome, &plus);
        let nt_m = crate::start_model::train_nterm_model(&rc, &minus);
        if let Some(ref model) = nt_p {
            crate::start_model::score_nterm_model(genome, &mut plus, model);
        }
        if let Some(ref model) = nt_m {
            crate::start_model::score_nterm_model(&rc, &mut minus, model);
        }
        let stop_p = crate::start_model::train_stop_model(genome, &plus);
        let stop_m = crate::start_model::train_stop_model(&rc, &minus);
        if let Some(ref model) = stop_p {
            crate::start_model::score_stop_model(genome, &mut plus, model);
        }
        if let Some(ref model) = stop_m {
            crate::start_model::score_stop_model(&rc, &mut minus, model);
        }
        for orf in plus.iter_mut() { orf.score = composite_score(orf, gc3_target); }
        for orf in minus.iter_mut() { orf.score = composite_score(orf, gc3_target); }
    }

    // 10. PWM-based RBS scoring
    let pwm_p = build_rbs_pwm(genome, &plus);
    let pwm_m = build_rbs_pwm(&rc, &minus);
    for orf in plus.iter_mut() {
        orf.rbs_pwm = score_rbs_pwm(genome, orf.seq_start, &pwm_p);
        orf.score = composite_score(orf, gc3_target);
    }
    for orf in minus.iter_mut() {
        orf.rbs_pwm = score_rbs_pwm(&rc, orf.seq_start, &pwm_m);
        orf.score = composite_score(orf, gc3_target);
    }

    // 10b. Compute upstream + downstream gene distance for start selection
    // Both distances inform operon context: real starts create normal spacing both ways.
    {
        // Plus strand: collect confident gene boundaries
        let mut plus_ends: Vec<usize> = plus.iter()
            .filter(|o| o.score > 0.45 && o.length >= 200 && o.is_longest)
            .map(|o| o.end).collect();
        let mut plus_starts: Vec<usize> = plus.iter()
            .filter(|o| o.score > 0.45 && o.length >= 200 && o.is_longest)
            .map(|o| o.start).collect();
        plus_ends.sort();
        plus_starts.sort();

        // Minus strand
        let mut minus_starts: Vec<usize> = minus.iter()
            .filter(|o| o.score > 0.45 && o.length >= 200 && o.is_longest)
            .map(|o| o.start).collect();
        let mut minus_ends: Vec<usize> = minus.iter()
            .filter(|o| o.score > 0.45 && o.length >= 200 && o.is_longest)
            .map(|o| o.end).collect();
        minus_starts.sort();
        minus_ends.sort();

        for orf in plus.iter_mut() {
            // Upstream: nearest gene end < our start
            let idx = plus_ends.partition_point(|&e| e < orf.start);
            if idx > 0 {
                orf.upstream_gene_dist = orf.start as i64 - plus_ends[idx - 1] as i64;
            }
            // Downstream: nearest gene start > our end (stored in stop_nn as temp, then moved)
            // Actually reuse upstream_gene_dist as min(upstream, downstream)
            let idx2 = plus_starts.partition_point(|&s| s <= orf.end);
            if idx2 < plus_starts.len() {
                let down_dist = plus_starts[idx2] as i64 - orf.end as i64;
                // Take minimum of upstream and downstream (closest neighbor)
                if orf.upstream_gene_dist == i64::MAX || down_dist < orf.upstream_gene_dist {
                    orf.upstream_gene_dist = down_dist;
                }
            }
        }

        for orf in minus.iter_mut() {
            // Minus strand: bio upstream = genomic downstream (higher coords)
            let idx = minus_starts.partition_point(|&s| s <= orf.end);
            if idx < minus_starts.len() {
                orf.upstream_gene_dist = minus_starts[idx] as i64 - orf.end as i64;
            }
            // Bio downstream = genomic upstream (lower coords)
            let idx2 = minus_ends.partition_point(|&e| e < orf.start);
            if idx2 > 0 {
                let down_dist = orf.start as i64 - minus_ends[idx2 - 1] as i64;
                if orf.upstream_gene_dist == i64::MAX || down_dist < orf.upstream_gene_dist {
                    orf.upstream_gene_dist = down_dist;
                }
            }
        }
    }

    // 11. Shadow detection on all ORFs
    {
        let mut all_shadow: Vec<Gene> = plus.iter().chain(minus.iter()).cloned().collect();
        detect_shadows(&mut all_shadow);
        let np = plus.len();
        for (i, orf) in all_shadow.iter().enumerate() {
            if i < np {
                plus[i].shadow_pen = orf.shadow_pen;
            } else {
                minus[i - np].shadow_pen = orf.shadow_pen;
            }
        }
    }

    // 12. Merge all ORFs for prediction
    let mut all_orfs: Vec<Gene> = Vec::new();
    all_orfs.extend(plus.iter().cloned());
    all_orfs.extend(minus.iter().cloned());

    let (mut results, _filtered) = run_prediction(&mut all_orfs, 0.0);

    // 13. Density adjustment (matching monolith)
    let total_coding: usize = results.iter().map(|&i| all_orfs[i].end - all_orfs[i].start + 1).sum();
    let density = total_coding as f64 / glen as f64;
    let density_adj: f64 = if density < 0.80 {
        -0.03
    } else if density > 0.93 {
        0.03
    } else {
        0.0
    };
    if density_adj.abs() > 0.001 {
        // Re-create fresh all_orfs (undo operon boosts from first run)
        let mut all_orfs2: Vec<Gene> = Vec::new();
        all_orfs2.extend(plus.iter().cloned());
        all_orfs2.extend(minus.iter().cloned());
        let (r2, _) = run_prediction(&mut all_orfs2, density_adj);
        results = r2;
        all_orfs = all_orfs2;
    }

    // 14. Second intergenic refinement (matching monolith)
    if results.len() > 200 {
        let res_orfs: Vec<Gene> = results.iter().map(|&i| all_orfs[i].clone()).collect();
        let ig2 = get_intergenic_regions(genome, &res_orfs, 50);
        if !ig2.is_empty() {
            let conf_p: Vec<Gene> = plus.iter()
                .filter(|o| o.score > 0.44 && o.length >= 250 && o.is_longest)
                .cloned().collect();
            let conf_m: Vec<Gene> = minus.iter()
                .filter(|o| o.score > 0.44 && o.length >= 250 && o.is_longest)
                .cloned().collect();

            let it2_p = train_intergenic_hex(genome, &conf_p, &ig2);
            let ig2_rc: Vec<Vec<u8>> = ig2.iter().map(|s| rev_comp(s)).collect();
            let it2_m = train_intergenic_hex(&rc, &conf_m, &ig2_rc);
            let it2 = merge_hex(&it2_p, &it2_m);

            if let Some(igm2) = it2 {
                hex_model = blend_hex(&hex_model, &igm2, 0.25);
                score_all_orfs(&mut plus, &mut minus, genome, &rc, &hex_model, &mono_model, gc3_target, true);

                // Re-build PWM
                let pwm_p2 = build_rbs_pwm(genome, &plus);
                let pwm_m2 = build_rbs_pwm(&rc, &minus);
                if pwm_p2.is_some() {
                    for orf in plus.iter_mut() {
                        orf.rbs_pwm = score_rbs_pwm(genome, orf.seq_start, &pwm_p2);
                        orf.score = composite_score(orf, gc3_target);
                    }
                }
                if pwm_m2.is_some() {
                    for orf in minus.iter_mut() {
                        orf.rbs_pwm = score_rbs_pwm(&rc, orf.seq_start, &pwm_m2);
                        orf.score = composite_score(orf, gc3_target);
                    }
                }

                // Re-detect shadows
                {
                    let mut all_s: Vec<Gene> = plus.iter().chain(minus.iter()).cloned().collect();
                    detect_shadows(&mut all_s);
                    let np = plus.len();
                    for (i, orf) in all_s.iter().enumerate() {
                        if i < np {
                            plus[i].shadow_pen = orf.shadow_pen;
                        } else {
                            minus[i - np].shadow_pen = orf.shadow_pen;
                        }
                    }
                }

                let mut all_orfs2: Vec<Gene> = Vec::new();
                all_orfs2.extend(plus.iter().cloned());
                all_orfs2.extend(minus.iter().cloned());
                let (results2, _) = run_prediction(&mut all_orfs2, density_adj);

                if results2.len() >= (results.len() as f64 * 0.85) as usize {
                    results = results2;
                    all_orfs = all_orfs2;
                }
            }
        }
    }

    // 14b. Iterative PWM retraining — DISABLED: neutral (F1=0.938 same as without)
    // Second-pass PWM trained on selected genes doesn't improve over first-pass.
    // Our first-pass training set (score>0.46, len>=350, ATG, longest) is already good.
    // Kept code + build_rbs_pwm_from_genes() in rbs.rs for future use.
    //
    // Original idea from Prodigal: retrain on DP results for better start selection.
    // May help on genomes with unusual SD usage (non-E.coli).
    // Use results from first DP as training set for second PWM pass.
    // Selected genes are higher quality than initial confident ORFs.
    if false && results.len() > 200 {  // DISABLED — see comment above
        let selected_genes: Vec<&Gene> = results.iter()
            .map(|&i| &all_orfs[i])
            .filter(|g| g.length >= 200 && g.is_atg())
            .collect();

        let pwm_p2 = crate::rbs::build_rbs_pwm_from_genes(genome, &selected_genes
            .iter().filter(|g| g.is_plus).copied().collect::<Vec<_>>(), 25);
        let pwm_m2 = crate::rbs::build_rbs_pwm_from_genes(&rc, &selected_genes
            .iter().filter(|g| !g.is_plus).copied().collect::<Vec<_>>(), 25);

        let mut rescored = false;
        if let Some(ref pwm) = pwm_p2 {
            for orf in plus.iter_mut() {
                let new_score = crate::rbs::score_rbs_pwm(genome, orf.seq_start, &Some(pwm.clone()));
                // Blend: 60% new PWM + 40% old
                // Use better of old vs new PWM score (not blend)
                orf.rbs_pwm = orf.rbs_pwm.max(new_score);
            }
            rescored = true;
        }
        if let Some(ref pwm) = pwm_m2 {
            for orf in minus.iter_mut() {
                let new_score = crate::rbs::score_rbs_pwm(&rc, orf.seq_start, &Some(pwm.clone()));
                // Use better of old vs new PWM score (not blend)
                orf.rbs_pwm = orf.rbs_pwm.max(new_score);
            }
            rescored = true;
        }
        if rescored {
            // Recompute composite scores
            for orf in plus.iter_mut() {
                orf.score = composite_score(orf, gc3_target);
            }
            for orf in minus.iter_mut() {
                orf.score = composite_score(orf, gc3_target);
            }
            // Re-run prediction with updated scores
            let mut all_orfs3: Vec<Gene> = Vec::new();
            all_orfs3.extend(plus.iter().cloned());
            all_orfs3.extend(minus.iter().cloned());
            let (r3, _) = run_prediction(&mut all_orfs3, density_adj);
            // Only accept if not much worse
            if r3.len() >= (results.len() as f64 * 0.90) as usize {
                results = r3;
                all_orfs = all_orfs3;
            }
        }
    }

    // 15. Atypical gene rescue
    rescue_atypical(&mut results, &all_orfs);

    // 15c. Targeted gap rescue (cycle 20): fill large intergenic gaps
    gap_targeted_rescue(&mut results, &all_orfs, glen, 400);

    // 15e. Operon rescue: genes inside operons don't need own RBS
    operon_rescue(&mut results, &all_orfs);

    // 15d. Filter same-strand overlaps > 45bp (cycle 20 cleanup, relaxed threshold)
    filter_same_strand_overlaps(&mut results, &all_orfs, 45);

    // 15f. Operon-internal gene rescue
    // If two same-strand selected genes have a gap >200bp, look for ORFs in that gap
    // with relaxed thresholds. Operon-internal genes don't need strong RBS/promoter.
    {
        let mut sel_sorted: Vec<(usize, usize, usize, bool)> = results.iter()
            .map(|&i| (all_orfs[i].start, all_orfs[i].end, i, all_orfs[i].is_plus))
            .collect();
        sel_sorted.sort();

        let sel_set: std::collections::HashSet<(usize, usize, bool)> = results.iter()
            .map(|&i| (all_orfs[i].start, all_orfs[i].end, all_orfs[i].is_plus)).collect();

        let mut rescued = 0;
        for w in sel_sorted.windows(2) {
            let (_, end1, _, strand1) = w[0];
            let (start2, _, _, strand2) = w[1];
            if strand1 != strand2 { continue; }
            let gap = if start2 > end1 { start2 - end1 } else { continue };
            if gap < 200 || gap > 2000 { continue; }

            // Find the BEST ORF in this gap (highest score, max 1 per gap)
            let mut best_idx: Option<usize> = None;
            let mut best_score = 0.0f64;

            for (i, orf) in all_orfs.iter().enumerate() {
                if sel_set.contains(&(orf.start, orf.end, orf.is_plus)) { continue; }
                if orf.is_plus != strand1 { continue; }
                if orf.start < end1.saturating_sub(10) { continue; }
                if orf.end > start2 + 10 { continue; }
                if orf.length < 100 { continue; }
                if !orf.is_longest { continue; }
                if !orf.is_atg() { continue; }
                if orf.score < 0.33 { continue; }
                // Need at least ONE strong positive signal
                if orf.hex_avg <= 0.05 && orf.rbs < 0.30 { continue; }

                if orf.score > best_score {
                    // Check no overlap
                    let mut ok = true;
                    for &ri in &results {
                        let a = &all_orfs[ri];
                        if orf.start.max(a.start) < orf.end.min(a.end) { ok = false; break; }
                    }
                    if ok {
                        best_score = orf.score;
                        best_idx = Some(i);
                    }
                }
            }

            if let Some(idx) = best_idx {
                results.push(idx);
                rescued += 1;
            }
        }
        if rescued > 0 {
            eprintln!("Operon-internal rescue: {} genes in same-strand gaps", rescued);
        }
    }

    // 15g. Leader peptide detection — DISABLED (ratio 7:204 TP:FP)
    // Domain detector works: finds thrL, ilvL, hisL, pheL, mgtL (5/12).
    // Three signals combined: AA enrichment + attenuator + operon context.
    // But too many FP short ORFs with palindromes near operon starts.
    // Needs: proper terminator detection (free energy) or Pfam TA families.
    // Code in start_model.rs: find_leader_peptides, detect_rho_independent_terminator.

    // 15g. Short gene rescue — DISABLED
    // Tested approaches:
    //   start_nn only (≥0.70, gaps): 326 rescued, ratio 18:288 TP:FP
    //   start_nn (≥0.90) + RBS: 76 rescued, ratio 9:62
    //   start_nn + stop_nn combined: 7 rescued, 0 TP
    // Root cause: self-supervised training on own predictions can't find
    // patterns for genes we DON'T predict. Need external training signal
    // (verified starts from experiments or comparative genomics).
    // Models kept: start_nn helps start selection (+0.001 F1), stop_nn computed but unused.

    // 15b. Refine starts — DISABLED (reduces F1 from 0.927 to 0.833)
    // results = refine_starts(&results, &all_orfs);

    // 16. Collect and deduplicate
    let mut output: Vec<Gene> = results.iter().map(|&i| all_orfs[i].clone()).collect();
    output.sort_by_key(|g| g.start);
    output.dedup_by(|a, b| a.start == b.start && a.end == b.end && a.is_plus == b.is_plus);

    // 17. Pseudogene detection — DISABLED
    // Tested 3 methods: frameshift, readthrough, fragment clustering.
    // All methods kill more TP than FP (ratio 5:1 to 8:1).
    // Root cause: pseudogenes are remnants of real genes — they have genuine
    // coding signal, RBS patterns, and start codons. Without homology data
    // (database search), they're indistinguishable from real genes.
    // 124 of 221 FP are pseudogenes, but we can't remove them without
    // also removing real genes. Code kept for future use with homology info.

    (output, all_orfs)
}


fn get_intergenic_regions(genome: &[u8], genes: &[Gene], min_gap: usize) -> Vec<Vec<u8>> {
    if genes.is_empty() { return vec![]; }
    let glen = genome.len();
    let mut intervals: Vec<(usize, usize)> = genes.iter()
        .map(|g| (g.start.saturating_sub(1), g.end))
        .collect();
    intervals.sort();
    let mut merged: Vec<(usize, usize)> = vec![intervals[0]];
    for &(s, e) in &intervals[1..] {
        let last = merged.last_mut().unwrap();
        if s <= last.1 { last.1 = last.1.max(e); } else { merged.push((s, e)); }
    }

    let mut regions = Vec::new();
    if merged[0].0 > min_gap {
        regions.push(genome[0..merged[0].0].to_vec());
    }
    for i in 0..merged.len() - 1 {
        let gs = merged[i].1;
        let ge = merged[i + 1].0;
        if ge > gs && ge - gs >= min_gap {
            regions.push(genome[gs..ge.min(glen)].to_vec());
        }
    }
    if glen > merged.last().unwrap().1 + min_gap {
        regions.push(genome[merged.last().unwrap().1..].to_vec());
    }
    regions
}

/// Parse GFF3 file with ncRNA regions (rRNA, tRNA).
/// Supports barrnap GFF3 output and tRNAscan-SE tabular output.
fn parse_ncrna_gff(path: &str) -> Vec<(usize, usize)> {
    let mut regions = Vec::new();
    let content = match fs::read_to_string(path) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Warning: cannot read ncRNA file {}: {}", path, e);
            return regions;
        }
    };
    for line in content.lines() {
        if line.starts_with('#') || line.starts_with("Sequence") || line.starts_with("---") || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 5 {
            // GFF3 format: seqid source type start end ...
            if let (Ok(start), Ok(end)) = (fields[3].trim().parse::<usize>(), fields[4].trim().parse::<usize>()) {
                regions.push((start.min(end), start.max(end)));
            }
        } else {
            // tRNAscan-SE tabular: Name tRNA# Begin End Type ...
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.len() >= 4 {
                if let (Ok(start), Ok(end)) = (fields[2].trim().parse::<usize>(), fields[3].trim().parse::<usize>()) {
                    regions.push((start.min(end), start.max(end)));
                }
            }
        }
    }
    regions
}

/// Filter genes that overlap significantly with ncRNA regions.
fn filter_ncrna_overlaps(genes: &mut Vec<Gene>, ncrna: &[(usize, usize)]) {
    if ncrna.is_empty() { return; }
    let before = genes.len();
    genes.retain(|g| {
        let glen = (g.end - g.start + 1) as f64;
        for &(ns, ne) in ncrna {
            let ov_s = g.start.max(ns);
            let ov_e = g.end.min(ne);
            if ov_e > ov_s {
                let overlap = (ov_e - ov_s) as f64;
                // Remove if >30% of gene overlaps ncRNA
                if overlap / glen > 0.30 {
                    return false;
                }
            }
        }
        true
    });
    let removed = before - genes.len();
    if removed > 0 {
        eprintln!("ncRNA filter: removed {} CDS predictions overlapping rRNA/tRNA", removed);
    }
}

// ═══════════════════════════════════════════════════════════════
// Pseudogene detection methods
// ═══════════════════════════════════════════════════════════════

/// Method 1: Frameshift detection.
/// A pseudogene fragment often has another ORF in a different frame on the same strand
/// covering the same region. If we find such frame-shifted ORF pairs, both are suspect.
/// Returns indices of genes to flag as pseudogene candidates.
fn detect_frameshift_fragments(genes: &[Gene], all_orfs: &[Gene]) -> Vec<bool> {
    let mut is_pseudo = vec![false; genes.len()];
    if genes.len() < 2 { return is_pseudo; }

    for (i, g) in genes.iter().enumerate() {
        // Look for other predicted genes on same strand, different frame, overlapping region
        for (j, h) in genes.iter().enumerate() {
            if i == j { continue; }
            if g.is_plus != h.is_plus { continue; }
            if g.frame == h.frame { continue; }

            let ov_s = g.start.max(h.start);
            let ov_e = g.end.min(h.end);
            if ov_e <= ov_s { continue; }

            let overlap = (ov_e - ov_s) as f64;
            let g_len = (g.end - g.start + 1) as f64;
            let h_len = (h.end - h.start + 1) as f64;

            // Both genes significantly overlap in different frames → frameshift signature
            let g_frac = overlap / g_len;
            let h_frac = overlap / h_len;

            if g_frac > 0.40 && h_frac > 0.40 {
                // Both are fragments — flag the shorter one
                if g_len <= h_len {
                    is_pseudo[i] = true;
                } else {
                    is_pseudo[j] = true;
                }
            } else if g_frac > 0.60 {
                // g is mostly contained in h's region but different frame
                is_pseudo[i] = true;
            }
        }
    }
    is_pseudo
}

/// Method 2: Read-through detection (premature stop codon).
/// Detects pseudogenes by finding coding signal continuing past the stop codon
/// in the same frame, BUT only if that region isn't covered by another predicted gene.
/// This distinguishes pseudogene readthrough from normal operon packing.
fn detect_readthrough(genes: &[Gene], genome: &[u8], rc: &[u8], hex: &HexModel) -> Vec<bool> {
    let mut is_pseudo = vec![false; genes.len()];

    for (i, g) in genes.iter().enumerate() {
        let seq = if g.is_plus { genome } else { rc };
        let se = g.seq_end;

        // Check: is the region past the stop covered by another predicted gene?
        let readthrough_start = g.end + 1;
        let readthrough_end = g.end + 300;
        let mut covered_by_other = false;
        for (j, h) in genes.iter().enumerate() {
            if i == j { continue; }
            if h.is_plus != g.is_plus { continue; }
            // Does h cover the readthrough region?
            let ov_s = readthrough_start.max(h.start);
            let ov_e = readthrough_end.min(h.end);
            if ov_e > ov_s {
                let coverage = (ov_e - ov_s) as f64 / (readthrough_end - readthrough_start).max(1) as f64;
                if coverage > 0.5 {
                    covered_by_other = true;
                    break;
                }
            }
        }
        if covered_by_other { continue; } // coding signal explained by another gene

        // Scan past stop codon in same frame
        let mut pos = se;
        let mut hex_sum = 0.0f64;
        let mut hex_n = 0u32;
        let mut first_start_offset = usize::MAX; // distance to first start codon

        while pos + 6 <= seq.len() && pos < se + 1500 {
            let codon = &seq[pos..pos+3];
            if crate::orf::is_stop(codon) { break; }
            if crate::orf::is_start(codon) && first_start_offset == usize::MAX {
                first_start_offset = pos - se;
            }
            if let Some(idx) = crate::io::hex_enc(&seq[pos..pos+6]) {
                hex_sum += hex[idx];
                hex_n += 1;
            }
            pos += 3;
        }

        // Key check: if there's a start codon early in the readthrough,
        // this is a normal ORF (next gene in operon), NOT a pseudogene.
        // Pseudogene readthrough has NO start codon — the frame just continues.
        if first_start_offset < 60 { continue; } // start within 60bp = new gene

        // Strong coding continuation with no start codon → premature stop → pseudogene
        if hex_n >= 30 {
            let readthrough_avg = hex_sum / hex_n as f64;
            let readthrough_bp = hex_n as usize * 3;
            if readthrough_avg > 0.3 && readthrough_bp as f64 / g.length.max(1) as f64 > 0.20 {
                is_pseudo[i] = true;
            }
        }
    }
    is_pseudo
}

/// Method 3: Fragment clustering (frameshift fragments).
/// Detects the specific pattern: two genes on same strand, tiny gap (< 10bp),
/// in different frames — the hallmark of a single insertion/deletion frameshift.
/// Only flags the SHORTER fragment. Very conservative to avoid hitting normal operons.
fn detect_fragment_clusters(genes: &[Gene], genome: &[u8], rc: &[u8], hex: &HexModel) -> Vec<bool> {
    let mut is_pseudo = vec![false; genes.len()];
    if genes.len() < 2 { return is_pseudo; }

    for i in 0..genes.len() - 1 {
        let g1 = &genes[i];

        for j in (i + 1)..genes.len() {
            let g2 = &genes[j];
            if g2.start > g1.end + 30 { break; }
            if g1.is_plus != g2.is_plus { continue; }
            if g1.frame == g2.frame { continue; } // same frame = not frameshift

            // Gap between the two genes (genomic coordinates)
            let gap = if g2.start > g1.end { g2.start - g1.end } else { 0 };
            let overlap = if g1.end >= g2.start { g1.end - g2.start + 1 } else { 0 };

            // Frameshift signature: very small gap (0-10bp) or tiny overlap (1-4bp)
            if gap > 10 && overlap == 0 { continue; }
            if overlap > 10 { continue; } // too much overlap = different feature

            let g1_len = g1.end - g1.start + 1;
            let g2_len = g2.end - g2.start + 1;

            // Both must be fragments: neither should be very long independently
            if g1_len > 3000 || g2_len > 3000 { continue; }

            // Both must have coding signal (they ARE from a real gene, just broken)
            if g1.hex_avg < 0.0 || g2.hex_avg < 0.0 { continue; }

            // Key check: operon genes have proper RBS, pseudogene fragments don't.
            // BOTH genes should have weak RBS — if either has good RBS, it's a real operon.
            if g1.rbs > 0.25 || g2.rbs > 0.25 { continue; }
            // Also check: downstream gene should not have strong start context
            let downstream = if g2.start > g1.start { g2 } else { g1 };
            if downstream.start_ctx > 0.20 { continue; }

            // Flag the shorter fragment
            if g1_len <= g2_len {
                is_pseudo[i] = true;
            } else {
                is_pseudo[j] = true;
            }
        }
    }
    is_pseudo
}

/// Apply pseudogene detection. Returns number removed.
/// mode: 1=frameshift, 2=readthrough, 3=fragments, 7=all combined
fn filter_pseudogenes(
    genes: &mut Vec<Gene>,
    all_orfs: &[Gene],
    genome: &[u8],
    rc: &[u8],
    hex: &HexModel,
    mode: u8,
) -> usize {
    let n = genes.len();
    let mut flagged = vec![false; n];

    if mode & 1 != 0 {
        let fs = detect_frameshift_fragments(genes, all_orfs);
        for i in 0..n { flagged[i] |= fs[i]; }
    }
    if mode & 2 != 0 {
        let rt = detect_readthrough(genes, genome, rc, hex);
        for i in 0..n { flagged[i] |= rt[i]; }
    }
    if mode & 4 != 0 {
        let fc = detect_fragment_clusters(genes, genome, rc, hex);
        for i in 0..n { flagged[i] |= fc[i]; }
    }

    let removed = flagged.iter().filter(|&&f| f).count();
    let mut idx = 0;
    genes.retain(|_| { let keep = !flagged[idx]; idx += 1; keep });
    removed
}
