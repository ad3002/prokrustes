# Prokrustes Pipeline: Technical Specification

## Overview

Prokrustes is a de novo bacterial gene finder. It takes a FASTA genome as input and outputs GFF3 gene predictions. No external databases required — all models are self-trained on the input genome.

Core idea: **iterative self-training** of hexamer coding models, followed by dynamic programming gene selection with multiple rescue and filtering passes.

```
Input FASTA
  → Step 0: Genetic code detection
  → Step 1: ORF finding (both strands)
  → Step 2: GC3 target computation
  → Step 3: RBS + upstream features
  → Step 3b: GC frame bootstrap (high-GC genomes only)
  → Step 4: Initial hexamer model training
  → Step 5: Monocodon model training
  → Step 6: Initial scoring
  → Step 7: Iterative self-training (12 rounds)
  → Step 8: Intergenic refinement
  → Step 9: Final scoring (hex + edge + HMM + neural models + PWM)
  → Step 11: Shadow detection
  → Step 12: DP gene selection (filter → start ranking → operon boost → weights → DP)
  → Step 13: Density adjustment
  → Step 14: Second intergenic refinement
  → Step 14c: Prophage detection (scoring law)
  → Step 15: Rescue passes (atypical, gap, operon-internal)
  → Step 15d-e: Overlap filters (same-strand 45bp, opposite-strand 200bp/50%)
  → Step 16: Collect + dedup
  → Step 16a: smORF rescue (gaps only, 90-300bp)
  → Step 16b: Data-driven confidence filter
  → Output GFF3
```

---

## Step 0: Genetic Code Detection

**Purpose:** Distinguish standard code 11 (TAA/TAG/TGA = stop) from code 4 (Mycoplasma: TGA = Trp).

**Method:**
1. Find top 200 longest ORFs under both codes
2. Compare median ORF length: ratio = median_code4 / median_standard

**Decision:**
- Code 4 if: ratio > 1.20 AND genome GC < 38%
- Otherwise: code 11

**Rationale:** All known code-4 organisms (Mycoplasma, Spiroplasma) are low-GC. In high-GC genomes TGA is naturally rare (AT-rich codon), inflating the ratio falsely.

---

## Step 1: ORF Finding

For each of 3 reading frames on each strand:
- Find all stop codons (in detected genetic code)
- Between consecutive stops: enumerate all start codons (ATG, GTG, TTG)
- Create candidate gene for each start → stop pair

**Minimum:** 90 bp (MIN_ORF).
Each ORF tagged with: `stop_group` ID, `is_longest` (longest start in its stop group), `frame`, `start_codon`.

---

## Step 2: GC3 Target

Compute genome-wide GC at 3rd codon position from confident long ORFs (≥900bp, ATG, longest). Used to normalize GC3 deviations in scoring.

Fallback chain: ≥900bp → ≥600bp → default 0.5.

---

## Step 3: RBS and Upstream Features

For each ORF, compute:

| Feature | Method | Range |
|---------|--------|-------|
| `rbs` | Shine-Dalgarno motif scan (15 variants, positions -24 to -3) | 0.0–1.0 |
| `upstream_at` | AT content 40bp upstream of start | 0.0–1.0 |
| `leaderless` | Promoter -10/-35 box detection | 0.0–1.0 |
| `start_ctx` | Kozak-like nucleotide pattern at -3/-1/+4 | 0.0–0.38 |
| `gc3` | GC at 3rd codon position | 0.0–1.0 |
| `gc3_bias` | |GC3 - GC12| wobble position bias | 0.0–1.0 |

---

## Step 3b: GC Frame Bootstrap (if gc3_target > 0.55)

**Problem:** In high-GC genomes, stop codons are rare → many long spurious ORFs → hexamer training contaminated.

**Solution:** Zero-knowledge coding signal before any training:
1. For each genome position: which of 3 frames has highest GC in 120bp window?
2. Real genes: consistent frame bias (wobble position). Random ORFs: inconsistent.
3. Compute genome-wide bias weights from long ORFs.
4. Score each ORF: dot product of per-ORF frame distribution with genome bias.
5. **Initial DP:** select genes using only length significance + GC frame score (no hexamers).
6. Train hexamers on this clean selected set instead of all long ORFs.

**Gate:** Only activated when gc3_target > 0.55 (Pseudomonas, M. tuberculosis, etc.)

---

## Step 4: Initial Hexamer Training

**Standard path (gc3 ≤ 0.55):**
- Training set: ORFs ≥900bp, longest in stop group
- Coding hexamers: in-frame (every position in ORF)
- Noncoding hexamers: out-of-frame + reverse strand
- Model: log-odds ratio per hexamer (4096 entries)
- Fallback: ≥600bp if first fails

**High-GC path (gc3 > 0.55):**
- Training set from GC frame DP selection (Step 3b)
- Uses `train_hex_from_set()` with min_len=100

---

## Step 5: Monocodon Tables

Train codon (3-mer) frequency model from ≥900bp ORFs. Separate plus/minus strand models, merged.

---

## Step 6: Initial Scoring

Score all ORFs with hex + mono models (no edge). Composite score computed (see Scoring section below). Filter: retain ORFs with score > 0.05.

Build hex index cache (HexCache) for fast rescoring in subsequent iterations.

---

## Step 7: Iterative Self-Training (12 rounds)

Each iteration i (0..11):
1. **Confidence threshold:** 0.46 - i × 0.02 (relaxes from 0.46 to 0.24)
2. **Min training length:** 450 - i × 30 (relaxes from 450 to 180)
3. Select ORFs passing threshold + length + is_longest
4. Train new hexamer model from selected set
5. **Blend:** weight = 0.35 + i × 0.05 (increases from 0.35 to 0.70)
6. Every 3rd iteration: also retrain monocodon
7. Rescore all ORFs with blended model (uses HexCache — no sequence re-scanning)

**Convergence check:** hash of first 50 hex values. Stop if delta < 0.1.

---

## Step 8: First Intergenic Refinement

**Condition:** ≥200 confident genes (score > 0.48, ≥350bp, longest).

1. Extract intergenic regions (gaps ≥50bp between confident genes)
2. Train separate hexamer on intergenic sequences
3. Blend with current model (weight 0.30)
4. Run 3 extra rounds at thresholds 0.40, 0.38, 0.36

---

## Step 9: Final Scoring

Full rescore with edge coding, then additional features:

**9a. HMM Viterbi:** Run full Viterbi gene finder, compute overlap fraction with each ORF → `viterbi_frac`.

**9b. Atypical gene rescue:** For low-GC genomes, train separate hex on AT-rich ORFs.

**9c. Neural models:** Start context NN, N-terminal composition, stop context → `start_nn`, `stop_nn`.

**9x. Upstream coding potential:** Score upstream region in same reading frame → `upstream_coding` (truncation evidence).

**10. RBS PWM:** Train position weight matrix from confirmed starts (score > 0.46, ≥350bp, ATG, longest). Rescore all ORFs → `rbs_pwm`.

**10b. Gene distance:** Compute distance to nearest upstream/downstream gene.

---

## Step 11: Shadow Detection

For each ORF, check if a stronger opposite-strand gene overshadows it:
- If opposite gene overlaps >65% AND has hex_avg > this + 0.5 AND is 1.2x longer → shadow_pen = 0.15
- If overlaps >45% with weaker criteria → shadow_pen = 0.35
- Otherwise shadow_pen = 1.0 (no penalty)

Shadow penalty multiplies DP weight (not score).

---

## Step 12: DP Gene Selection (`run_prediction`)

### 12.1 Filter

Hard filters (rejected regardless of score):
```
hex_avg < -0.3 AND length < 600 → REJECT
length < 250 AND hex_avg < 0.05 AND rbs < 0.15 → REJECT
frame_bias < -0.4 AND length < 500 → REJECT
length < 180 AND fb < 0.10 AND hex < 0.20 AND rbs < 0.40 → REJECT
```

Score threshold (length-dependent):
```
length ≥ 900: thresh = 0.10
length ≥ 600: thresh = 0.19
length ≥ 300: thresh = 0.31
length ≥ 150: thresh = 0.39
else:         thresh = 0.49

thresh += density_adj           # ±0.03 from Step 13
thresh -= region_bonus          # prophage routing (Step 14c)
if rbs > 0.60: thresh *= 0.80  # strong RBS lowers bar
if rbs > 0.35: thresh *= 0.90

KEEP if score ≥ thresh
```

### 12.2 Keep Best 4 Starts per Stop Group

For each stop group, rank starts by:
```
sq = 0.18*rbs_combined + 0.09*hex_norm + 0.06*frame_bias + 0.06*edge
   + 0.14*start_type + 0.12*(1 - exp(-len/500))
   + nn_bonus + nterm_bonus + longest_bonus - truncation_penalty
```
Keep top 4.

### 12.3 Operon Boost

For each gene: sum scores of same-strand neighbors within 250bp.
- Sum > 1.5 → +0.08
- Sum > 0.8 → +0.05
- Sum > 0.4 → +0.03

### 12.4 Weight Computation

```
base = {0.12, 0.21, 0.33, 0.41, 0.51} by length bracket
base += density_adj
if rbs > 0.60: base *= 0.80
if rbs > 0.35: base *= 0.90

weight = (score - base) + max(0, hex_total * 0.004)
```

**Short gene penalty asymmetry** (for length < 250bp):
- If weight < 0: weight *= 250/length (amplify negatives)
- If weight > 0: weight *= length/250 (shrink positives)

### 12.4b Connection Score

Same-strand neighbor bonus added to weight:
- Overlap 1-4bp (ATGA coupling): +0.04
- Gap 0-20bp: +0.03
- Gap 20-50bp: +0.02
- Gap 50-150bp: +0.01

### 12.5-6 Dynamic Programming

Per-strand iterative DP (4 iterations, max_overlap=60bp), then combined DP (max_overlap=50bp). Shadow penalty applied as weight multiplier.

### 12.7 Gap Fill

Fill gaps with candidates: score ≥ 0.40, length ≥ 120, overlap ≤ 45bp (same strand) / 90bp (opposite).

---

## Step 13: Density Adjustment

```
density = total_coding_bp / genome_len
if density < 0.80: adj = -0.03  (relax thresholds)
if density > 0.93: adj = +0.03  (tighten)
```
Re-run prediction if adjusted.

---

## Step 14: Second Intergenic Refinement

Same as Step 8 but using DP results. Accept if ≥85% of genes retained.

---

## Step 14c: Prophage Detection

**Scoring law** (calibrated from 31 TP + 247 FP segments across 6 genomes):

1. Compute hex_threshold = 10th percentile of genes ≥600bp
2. Find clusters of low-hex genes (allowing 1-2 gaps)
3. For each candidate segment compute:
   - F = lowHexFraction / 0.5 (saturates at 50%)
   - D = density_excess (segment vs genome), positive only
   - size_score = span_kb / 20
   - gene_score = n_genes / 20
   - frag_penalty = fragmentation (0 = one block, 1 = all singletons)

```
Score = 0.8·F + 0.6·D + 0.2·size + 1.0·n_genes - 0.8·fragmentation
Threshold = 2.3
```

**Routing:** Currently only weak penalty (-0.03) in degraded regions. Prophage regions: no bonus (precision 64% insufficient for safe routing).

---

## Step 15: Rescue Passes

**15a. Atypical gene rescue:** Strong RBS + moderate coding.
**15c. Gap rescue:** Fill gaps 400-3000bp with confident candidates (score > 0.42).
**15e. Operon-internal rescue:** Fill same-strand gaps 200-2000bp (score > 0.33).
**15d. Same-strand overlap filter:** Remove if overlap > 45bp.
**15e. Opposite-strand overlap filter:** Remove if overlap > 200bp OR > 50% of either gene.

---

## Step 16: Output Assembly

**16. Collect + dedup** sorted by position.

**16a. smORF rescue:** In intergenic gaps, find short ORFs (90-300bp) with score > 0.33, hex > 0.05, and (RBS > 0.30 OR hex > 0.15).

**16b. Confidence filter** (data-driven, from 35K TP / 3.5K FP):

| Length | Precision | Required evidence |
|--------|-----------|-------------------|
| ≥600bp | 97% | Always keep |
| 300-600bp | 89% | hex > 0 OR rbs > 0.30 OR frame_bias > 0.10 |
| 150-300bp | 71% | hex > 0.05 OR (rbs > 0.35 AND is_longest) |
| <150bp | 18% | (hex > 0.10 AND rbs > 0.20) OR (hex > 0.15 AND longest) OR (rbs > 0.50 AND fb > 0.10) |

**Output:** GFF3 with seqid, prokrustes, CDS, start, end, score, strand, 0, attributes (ID, start_type, rbs_score, hex_score).

---

## Composite Score Function

```
score = 0.28 * hex_norm          // hexamer coding potential (dominant)
      + 0.15 * len_norm          // 1 - exp(-len/300)
      + 0.14 * rbs_combined      // max(rbs, 0.85 * rbs_pwm_norm)
      + 0.08 * longest_bonus     // 0.06 if is_longest
      + 0.06 * at_bonus          // upstream AT vs genome background
      + 0.05 * start_type        // ATG=1.0, GTG=0.55, TTG=0.35
      + 0.04 * frame_bias_norm   // codon frequency asymmetry
      + 0.04 * hex_cov           // fraction with positive hex
      + 0.03 * gc3_norm          // GC3 deviation from target
      + 0.01 * edge_norm         // upstream/downstream coding contrast
      + bonuses (RBS rescue +0.04, leaderless +0.05, start_ctx +0.02, gc3_bias +0.02)
```

Weights from logistic regression on 192K ORFs (E. coli training data).

---

## Disabled Features

| Feature | Tested F1 Impact | Reason disabled |
|---------|-----------------|-----------------|
| Leader peptide detection | 7:204 TP:FP | Needs terminator detection |
| Short gene rescue (NN) | 18:288 TP:FP | Self-training can't find absent genes |
| Start refinement | -0.094 F1 | Massive regression |
| Pseudogene detection | 5:1 kill ratio | Indistinguishable without homology |
| Iterative PWM retraining | 0.000 | First-pass already optimal |
| Dicodon scoring | 0.000 | Correlated with hex, no new info |
| Upstream composition model | 0.000 | Too small effect at safe coefficients |

---

## Multi-Contig Support

For multi-record FASTA (e.g. Borrelia with 22 replicons):
- Each contig annotated independently (own hex training, own DP)
- GFF output uses correct seqid per contig
- Short contigs (<300bp) skipped

---

## Performance

10-genome benchmark (with ncRNA masking):

| Genome | GC% | F1 | vs Prodigal |
|--------|-----|-----|-------------|
| E. coli K-12 | 51 | 0.946 | -0.004 |
| Salmonella LT2 | 52 | 0.951 | -0.005 |
| Bacillus subtilis | 44 | 0.962 | -0.006 |
| S. aureus NCTC 8325 | 33 | 0.943 | -0.004 |
| **Mycoplasma genitalium** | 32 | **0.937** | **+0.689** |
| Borrelia burgdorferi | 29 | 0.918 | -0.013 |
| Bacteroides fragilis | 43 | 0.944 | -0.015 |
| Synechocystis PCC 6803 | 48 | 0.931 | -0.023 |
| Pseudomonas aeruginosa | 67 | 0.928 | -0.049 |
| M. tuberculosis H37Rv | 66 | 0.876 | -0.052 |
