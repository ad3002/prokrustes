#!/usr/bin/env bash
set -euo pipefail

# Prokrustes — one-command setup
# Builds the binary, downloads reference data, installs ncRNA tools
#
# Usage:
#   ./scripts/setup.sh           # full setup
#   ./scripts/setup.sh --no-bio  # skip barrnap/tRNAscan (Rust only)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$PROJECT_DIR/data"

echo "=== Prokrustes Setup ==="
echo ""

# --- 1. Build ---
echo "[1/4] Building prokrustes..."
cd "$PROJECT_DIR"
if ! command -v cargo &>/dev/null; then
    echo "  Rust not found. Installing via rustup..."
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
    source "$HOME/.cargo/env"
fi
cargo build --release 2>&1 | tail -3
echo "  -> target/release/prokrustes built"

# --- 2. Download reference genomes ---
echo ""
echo "[2/4] Downloading reference genomes..."
mkdir -p "$DATA_DIR"

download() {
    local url="$1" dest="$2" label="$3"
    if [ -f "$dest" ]; then
        echo "  $label — already exists"
    else
        echo "  $label — downloading..."
        curl -sL "$url" | gunzip > "$dest"
        echo "  -> $(wc -c < "$dest" | tr -d ' ') bytes"
    fi
}

download \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" \
    "$DATA_DIR/ecoli_k12.fasta" \
    "E. coli K-12 MG1655 genome"

download \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz" \
    "$DATA_DIR/ecoli_k12.gff" \
    "E. coli K-12 reference annotations"

download \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz" \
    "$DATA_DIR/salmonella_lt2.fasta" \
    "Salmonella LT2 genome"

download \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz" \
    "$DATA_DIR/salmonella_lt2.gff" \
    "Salmonella LT2 reference annotations"

download \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz" \
    "$DATA_DIR/bacillus_subtilis.fasta" \
    "Bacillus subtilis 168 genome"

# --- 3. Install bioinformatics tools (optional) ---
if [[ "${1:-}" != "--no-bio" ]]; then
    echo ""
    echo "[3/4] Installing ncRNA annotation tools..."

    if command -v conda &>/dev/null || command -v mamba &>/dev/null; then
        CONDA_CMD=$(command -v mamba || command -v conda)
        if ! command -v barrnap &>/dev/null; then
            echo "  Installing barrnap via $CONDA_CMD..."
            $CONDA_CMD install -y -c bioconda barrnap 2>&1 | tail -1
        else
            echo "  barrnap — already installed"
        fi
        if ! command -v tRNAscan-SE &>/dev/null; then
            echo "  Installing tRNAscan-SE via $CONDA_CMD..."
            $CONDA_CMD install -y -c bioconda trnascan-se 2>&1 | tail -1
        else
            echo "  tRNAscan-SE — already installed"
        fi
    else
        echo "  conda/mamba not found — skipping ncRNA tools"
        echo "  Install manually: conda install -c bioconda barrnap trnascan-se"
    fi

    # Generate ncRNA annotations if tools available
    if command -v barrnap &>/dev/null; then
        if [ ! -f "$DATA_DIR/ecoli_k12_ncrna.gff" ]; then
            echo "  Generating E. coli ncRNA annotations..."
            barrnap "$DATA_DIR/ecoli_k12.fasta" > "$DATA_DIR/ecoli_k12_ncrna.gff" 2>/dev/null
        fi
    fi
else
    echo ""
    echo "[3/4] Skipping ncRNA tools (--no-bio)"
fi

# --- 4. Verify ---
echo ""
echo "[4/4] Verification..."
BINARY="$PROJECT_DIR/target/release/prokrustes"

echo "  Running quick annotation test..."
PRED=$("$BINARY" "$DATA_DIR/ecoli_k12.fasta" 2>/dev/null | grep -c $'\tCDS\t' || echo 0)
echo "  -> Predicted $PRED CDS on E. coli K-12"

if [ "$PRED" -gt 3000 ]; then
    echo ""
    echo "=== Setup complete ==="
    echo ""
    echo "Quick start:"
    echo "  # Annotate a genome"
    echo "  target/release/prokrustes data/ecoli_k12.fasta > output.gff"
    echo ""
    echo "  # With ncRNA masking"
    echo "  target/release/prokrustes data/ecoli_k12.fasta --ncrna data/ecoli_k12_ncrna.gff > output.gff"
    echo ""
    echo "  # With comparative genomics"
    echo "  target/release/prokrustes data/ecoli_k12.fasta --compare data/salmonella_lt2.fasta > output.gff"
    echo ""
    echo "  # Run benchmark"
    echo "  bash scripts/benchmark.sh"
else
    echo ""
    echo "WARNING: only $PRED CDS predicted — something may be wrong"
    exit 1
fi
