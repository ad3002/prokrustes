# Prokrustes

A fast prokaryotic genome annotation tool built in Rust. Prokrustes uses epistemic layering — routing genomic regions to specialized scoring contexts — to achieve near-Prodigal accuracy with fully interpretable error attribution.

## Features

- **Three-layer architecture**: region classification, transcription unit context, per-gene interpretation
- **LightGBM start codon ranking**: learning-to-rank on 22 features (RBS, Kozak-like context, upstream AT, GC3 bias, etc.)
- **IS element detection**: terminal inverted repeat (TIR) + direct repeat scanning
- **Prophage-aware mode**: alien codon priors for horizontally transferred regions
- **Intrinsic terminator finder**: Rho-independent hairpin detection with thermodynamic scoring
- **Conservation scoring**: Dayhoff-6 reduced alphabet k-mer comparison across related genomes
- **HMM gene finding**: hexanucleotide Viterbi with self-training (3-pass blend)
- **Zero external runtime dependencies**: pure Rust inference, no ML frameworks needed at runtime

## Performance

Primary benchmark: *E. coli* K-12 MG1655 (4.6 Mb, ~4400 CDS)

| Configuration | Precision | Recall | F1 |
|---|---|---|---|
| Prodigal (baseline) | 0.957 | 0.943 | **0.950** |
| Prokrustes (current) | 0.951 | 0.944 | **0.948** |

Gap: 0.002 F1. Validated on 10 prokaryotic genomes across diverse clades.

## Quick Start

### Option 1: One-command setup (recommended)

```bash
git clone https://github.com/ad3002/prokrustes.git
cd prokrustes
./scripts/setup.sh
```

This will:
- Build the Rust binary
- Download 10 reference genomes from NCBI
- Install ncRNA tools (barrnap, tRNAscan-SE) if conda is available
- Run a verification test on E. coli K-12

### Option 2: Docker (fully reproducible, no dependencies)

```bash
docker build -t prokrustes .

# Full benchmark: 10 genomes, Prokrustes vs Prodigal, prints comparison table
docker run prokrustes benchmark

# Annotate a single genome
docker run prokrustes annotate /data/ecoli_k12.fasta > output.gff

# See all options
docker run prokrustes --help
```

The Docker image includes all 10 reference genomes, Prodigal baseline, and ncRNA tools.

### Option 3: Manual build

Requirements: Rust 1.70+

```bash
git clone https://github.com/ad3002/prokrustes.git
cd prokrustes
cargo build --release
```

The binary will be at `target/release/prokrustes`.

## Reproducible Benchmark

Run the full benchmark across 10 genomes with Prodigal comparison:

```bash
# Install Prodigal baseline (optional)
conda install -c bioconda prodigal
# or: apt install prodigal

# Run full benchmark
./scripts/benchmark.sh
```

This downloads all genomes, runs both Prokrustes and Prodigal, and prints a comparison table:

```
Genome                       Size   Prokrustes F1    Prodigal F1    Delta
---------------------------------------------------------------------------
Escherichia coli K-12         4.6          0.9480         0.9500  -0.0020
Salmonella enterica LT2       5.0            ...            ...      ...
Bacillus subtilis 168          4.2            ...            ...      ...
...
```

Run a single genome:

```bash
./scripts/benchmark.sh ecoli_k12
```

## Usage

### Basic annotation

```bash
prokrustes genome.fasta
```

Output: GFF3-format gene predictions to stdout, statistics to stderr.

### With ncRNA masking

To exclude rRNA/tRNA regions from gene prediction, provide a GFF file (e.g., from [barrnap](https://github.com/tseemann/barrnap) or [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/)):

```bash
# Generate ncRNA annotations
barrnap genome.fasta > ncrna.gff
tRNAscan-SE -B genome.fasta >> ncrna.gff

# Run with masking
prokrustes genome.fasta --ncrna ncrna.gff
```

### With comparative genomics

Provide one or more related genomes to improve predictions via conservation scoring:

```bash
prokrustes genome.fasta --compare related1.fasta related2.fasta
```

### Debug a specific region

```bash
prokrustes genome.fasta --debug <start> <end> <strand>
# Example: prokrustes genome.fasta --debug 100000 105000 + 2
```

### Additional output modes

```bash
# Dump all candidate ORFs (not just selected genes)
prokrustes genome.fasta --dump-all-orfs

# Export feature matrix as TSV
prokrustes genome.fasta --dump-tsv

# HMM-only mode (no start ranker)
prokrustes genome.fasta --hmm
```

## Input Data

### Genome sequences

Download reference genomes from NCBI:

```bash
# E. coli K-12 MG1655 (primary benchmark)
curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" | gunzip > ecoli_k12.fasta

# Reference annotations (for benchmarking)
curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz" | gunzip > ecoli_k12.gff

# Salmonella LT2 (validation)
curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz" | gunzip > salmonella_lt2.fasta

# Bacillus subtilis 168 (validation)
curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz" | gunzip > bacillus_subtilis.fasta
```

### ncRNA tools (optional, recommended)

```bash
# Install via conda
conda install -c bioconda barrnap trnascan-se
```

### LightGBM model

The pre-trained start codon ranker is included at `models/start_ranker.lgb`. It was trained on *E. coli* K-12 with 4146 loci and 84907 candidate starts.

To retrain on your own data, see `python/prokrustes/pipeline.py`.

## Architecture

```
Layer 1 — Region Classification (5-10 Kb windows)
  ├── Core bacterial genome → standard pipeline
  ├── IS/mobile elements → TIR + repeat detection, relaxed overlap rules
  ├── Prophage islands → alien-codon priors, permissive thresholds
  └── Intergenic deserts → smORF-specific search mode

Layer 2 — Transcription Unit Context
  ├── Promoter → terminator boundaries (operon architecture)
  ├── Gene position: first / internal / last in operon
  ├── Overlap & coupling patterns
  └── Conservation signal (Dayhoff-6 reduced alphabet, multi-genome panel)

Layer 3 — Per-Gene Interpretation
  ├── Main CDS detection (hexamer HMM + Viterbi)
  ├── Start selection (LightGBM lambdarank on 22 features)
  ├── smORF pipeline (separate short-gene ranker)
  ├── Leader peptide detection (operon context + AA enrichment)
  └── Pseudogene flagging (conservation + structural breakage)
```

## Output Format

Standard GFF3 to stdout:

```
##gff-version 3
seqid  prokrustes  CDS  start  end  score  strand  phase  ID=gene_N;rbs_score=X;start_type=ATG;hex_score=Y
```

## Known Limitations

- **Short genes (<150 bp)**: hexamer signal too noisy; smORF pipeline in development
- **Pseudogenes**: cannot distinguish from real genes without homology data
- **Mycoplasma / low-GC genomes**: requires multi-clade training (known gap)
- **Start accuracy**: ~30% of RefSeq starts may be incorrectly annotated, complicating validation

## Python Training Bridge

The `python/prokrustes/` directory contains utilities for:
- Training the LightGBM start ranker (`pipeline.py`)
- Running the Rust binary and parsing output (`pipeline_rust.py`)
- Evaluating predictions against reference annotations

```bash
pip install lightgbm pandas scikit-learn
python -m prokrustes.pipeline --train --genome ecoli_k12.fasta --reference ecoli_k12.gff
```

## Citation

If you use Prokrustes in your research, please cite:

> Prokrustes: prokaryotic genome annotation via epistemic layering. (2026). https://github.com/ad3002/prokrustes

## License

MIT
