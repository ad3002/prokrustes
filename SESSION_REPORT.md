# Prokrustes Session Report — April 5-6, 2026

## Gene Finding Performance

| Genome | Before | After | Prodigal | vs Prodigal |
|--------|--------|-------|----------|-------------|
| **Mycoplasma genitalium** | 0.925 | **0.949** | 0.248 | **WIN +0.700** |
| Bacillus subtilis | 0.958 | **0.965** | 0.968 | -0.004 |
| Salmonella LT2 | 0.949 | **0.950** | 0.956 | -0.006 |
| E. coli K-12 | 0.943 | **0.945** | 0.950 | -0.005 |
| S. aureus NCTC 8325 | 0.943 | **0.944** | 0.947 | -0.003 |
| Bacteroides fragilis | 0.927 | **0.949** | 0.959 | -0.010 |
| Synechocystis PCC 6803 | 0.831 | **0.936** | 0.954 | -0.018 |
| **Pseudomonas aeruginosa** | **0.230** | **0.945** | 0.976 | -0.032 |
| Borrelia burgdorferi | 0.537 | **0.921** | 0.931 | -0.011 |
| **M. tuberculosis H37Rv** | **0.443** | **0.884** | 0.928 | -0.044 |

## What Worked (by impact)

### 1. V2 Start Ranker (+0.005 to +0.017 F1) — LARGEST SINGLE IMPROVEMENT
- Trained LightGBM pairwise ranker from 26006 stop groups (537K samples, 22 features)
- NDCG@1 = 0.888 on validation, top-1 accuracy 88%
- Data showed: model over-preferred is_longest (47% vs 26% correct) and RBS (64% wrong had better RBS)
- V2 model learned correct feature balance from data instead of hand-tuned weights
- Top features: frac_of_longest, length_rank, delta_score_to_best, start_type

### 2. GC Frame Bootstrap (+0.71 on Pseudomonas)
- High-GC genomes (>55% GC): hexamer training contaminated by spurious long ORFs
- Solution: zero-knowledge GC frame bias → initial DP without hexamers → clean training set
- Also fixed genetic code detection (GC guard prevents false code-4 on high-GC)

### 3. Multi-contig FASTA (+0.38 on Borrelia)
- Borrelia has 22 replicons — concatenating them broke coordinates for all genes after first contig
- Each contig now annotated independently with correct seqid in GFF output

### 4. Data-driven Confidence Filter
- Learned thresholds from 35K TP / 3.5K FP across 10 genomes
- Key: <150bp genes have 18% precision → strict multi-signal requirement
- hex_avg is best single discriminator (TP 5th pct = 0.141, FP median = 0.107)

### 5. Structural Filters
- Opposite-strand overlap limit (200bp / 50% — Prodigal standard)
- Short gene penalty asymmetry (amplify negatives, shrink positives for <250bp)
- ncRNA masking (removes ~30 FP per genome from rRNA/tRNA regions)

## What Didn't Work

### Score Formula Changes
- Log-odds length factor, GC-adjusted len_scale, upstream composition scoring
- Any change to composite_score disrupts the weight balance trained on E. coli data
- Lesson: weights were ML-calibrated together; changing one breaks all

### Prophage Routing
- Detection precision maxed at 64% (from 8%) with calibrated scoring law
- Dual hex model: +0.017 phage-region F1, but overall neutral due to FP regions
- Integrase markers: too diverse for de novo detection (81 sequences → 68 singleton families)
- att repeat detection: candidate boundaries too imprecise (±10-20kb vs real ±100bp)
- **Conclusion**: needs external DB (Pfam integrases) for >90% precision

### Dicodon Scoring
- Correlated with hexamer — same information from same data, no improvement

### Expression-level Dual Hex
- 100% of genome blocks have different codon usage (expression-level bias)
- But dual hex model was ~neutral: V2 ranker already handles the variation

## Synteny Infrastructure (Built, Ready for Next Phase)

### Anchor Graph
- 832 anchor genes (present in ≥4/10 genomes)
- 3582 edges (anchor-to-anchor intervals)
- 854 conserved edges (≥2 genomes)

### Anchor Library
- 100 universal single-copy protein families
- 4-mer fingerprint profiles per family
- Leave-one-out validated: **96.6% precision, 69.1% recall**

### Genome Normalization
- All 10 genomes rotated to dnaA origin on (+) strand
- E. coli ↔ Salmonella: 87% same-strand synteny, 48 collinear blocks

### Analysis Results
- Anchor intervals: 85% have 0 expected genes (operonic, too dense for interval rescue)
- Expansion detector: 15 lineage-specific gene insertions (not phage)
- Local codon variation: 100% of blocks differ from genome average (expression bias)

### Interactive Visualization
- D3.js circular chromosome viewer at `viz/index.html`
- Zoom/pan, genome selection, node inspection, edge filtering

## Key Methodological Insights

1. **Data first, code second**: measuring phage error concentration (3x in 5% genome) was worth more than 20 iterations of threshold tuning
2. **ML > heuristics for start selection**: 1463 labeled pairs → LightGBM → 88% accuracy vs hand-tuned formula
3. **Calibrate from data, not from head**: prophage scoring law F1 jumped from ~0.20 to 0.65 with grid-searched weights
4. **Don't change pre-DP weights globally**: cascading effects through operon boost, gap fill, connection scoring are unpredictable
5. **Precision > recall for routing**: 64% precision makes any downstream intervention net-negative

## Architecture (PIPELINE.md)

Full 17-step pipeline documented in PIPELINE.md with every threshold, weight, condition, and disabled feature with rationale.

## Next Steps (Planned, Not Started)

1. **Anchor-chain placer**: insert new genome into graph by concordant anchor chains
2. **Anchor-guided annotation**: use expected gene count between anchors for rescue/pruning
3. **Local coding models**: train block-conditioned hex from anchor-defined neighborhoods
4. **Expansion-aware annotation**: special handling for lineage-specific insertion regions

## Repository

- GitHub: https://github.com/ad3002/prokrustes
- Forum post: https://forum.omniscale.space/t/prokaryotic-gene-finder-closing-the-last-0-6-to-prodigal-via-epistemic-layering/26
- Docker: `docker build -t prokrustes . && docker run prokrustes benchmark`
