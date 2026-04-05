#!/usr/bin/env bash
set -euo pipefail

# Prokrustes benchmark script
# Downloads E. coli K-12 reference, runs annotation, computes F1

DATA_DIR="${1:-data}"
BINARY="${2:-target/release/prokrustes}"

mkdir -p "$DATA_DIR"

# --- Download reference genome and annotations ---
GENOME="$DATA_DIR/ecoli_k12.fasta"
REFERENCE="$DATA_DIR/ecoli_k12.gff"

if [ ! -f "$GENOME" ]; then
    echo "Downloading E. coli K-12 MG1655 genome..."
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" | gunzip > "$GENOME"
    echo "  -> $(wc -c < "$GENOME") bytes"
fi

if [ ! -f "$REFERENCE" ]; then
    echo "Downloading E. coli K-12 MG1655 reference annotations..."
    curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz" | gunzip > "$REFERENCE"
    echo "  -> $(wc -c < "$REFERENCE") bytes"
fi

# --- Build if needed ---
if [ ! -f "$BINARY" ]; then
    echo "Building prokrustes (release)..."
    cargo build --release 2>&1 | tail -3
fi

# --- Run annotation ---
echo ""
echo "Running Prokrustes on E. coli K-12..."
PREDICTION="$DATA_DIR/prokrustes_prediction.gff"
time "$BINARY" "$GENOME" > "$PREDICTION" 2> "$DATA_DIR/prokrustes_stderr.log"

PRED_COUNT=$(grep -c $'\tCDS\t' "$PREDICTION" || echo 0)
echo "  Predicted: $PRED_COUNT CDS features"

# --- Evaluate ---
echo ""
echo "Computing F1 against NCBI reference..."
python3 "$(dirname "$0")/evaluate.py" "$PREDICTION" "$REFERENCE"
