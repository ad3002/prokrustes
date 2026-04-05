#!/usr/bin/env bash
set -euo pipefail

# Prokrustes — full reproducible benchmark with Prodigal baseline
# Downloads 10 prokaryotic genomes, annotates each with Prokrustes and Prodigal,
# evaluates F1 for both, prints comparison table.
#
# Usage:
#   ./scripts/benchmark.sh                  # all genomes
#   ./scripts/benchmark.sh ecoli_k12        # single genome
#   ./scripts/benchmark.sh --no-prodigal    # skip Prodigal baseline

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Auto-detect Docker vs local environment
if [ -f /usr/local/bin/prokrustes ] && [ -d /data ]; then
    # Docker
    DATA_DIR="/data"
    RESULTS_DIR="/results"
    BINARY="/usr/local/bin/prokrustes"
else
    # Local
    DATA_DIR="$PROJECT_DIR/data"
    RESULTS_DIR="$PROJECT_DIR/results"
    BINARY="$PROJECT_DIR/target/release/prokrustes"
fi

mkdir -p "$DATA_DIR" "$RESULTS_DIR"

# --- Genome registry ---
# Format: key|label|fasta_url|gff_url
GENOMES=(
    "ecoli_k12|Escherichia coli K-12 MG1655|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz"
    "salmonella_lt2|Salmonella enterica LT2|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz"
    "bacillus_subtilis|Bacillus subtilis 168|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.gff.gz"
    "pseudomonas_aeruginosa|Pseudomonas aeruginosa PAO1|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.gff.gz"
    "saureus_nctc8325|Staphylococcus aureus NCTC 8325|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz"
    "mycobacterium_tuberculosis|Mycobacterium tuberculosis H37Rv|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz"
    "mycoplasma_genitalium|Mycoplasma genitalium G37|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz"
    "synechocystis_pcc6803|Synechocystis sp. PCC 6803|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/725/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.fna.gz|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/725/GCF_000009725.1_ASM972v1/GCF_000009725.1_ASM972v1_genomic.gff.gz"
    "bacteroides_fragilis|Bacteroides fragilis YCH46|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/925/GCF_000009925.1_ASM992v1/GCF_000009925.1_ASM992v1_genomic.fna.gz|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/925/GCF_000009925.1_ASM992v1/GCF_000009925.1_ASM992v1_genomic.gff.gz"
    "borrelia_burgdorferi|Borrelia burgdorferi B31|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/685/GCF_000008685.2_ASM868v2/GCF_000008685.2_ASM868v2_genomic.fna.gz|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/685/GCF_000008685.2_ASM868v2/GCF_000008685.2_ASM868v2_genomic.gff.gz"
)

# --- Parse arguments ---
FILTER=""
RUN_PRODIGAL=true
for arg in "$@"; do
    if [ "$arg" = "--no-prodigal" ]; then
        RUN_PRODIGAL=false
    elif [[ "$arg" != --* ]]; then
        FILTER="$arg"
    fi
done

# --- Check Prodigal ---
PRODIGAL_BIN=""
if $RUN_PRODIGAL; then
    if command -v prodigal &>/dev/null; then
        PRODIGAL_BIN="prodigal"
    else
        echo "NOTE: Prodigal not found."
        echo "      Install: conda install -c bioconda prodigal"
        echo "      Or:      apt install prodigal"
        echo "      Running without baseline comparison."
        echo ""
        RUN_PRODIGAL=false
    fi
fi

# --- Build Prokrustes if needed ---
if [ ! -f "$BINARY" ]; then
    echo "Building prokrustes (release)..."
    cd "$PROJECT_DIR" && cargo build --release 2>&1 | tail -3
fi

# --- Download helper ---
download() {
    local url="$1" dest="$2"
    if [ ! -f "$dest" ]; then
        curl -sL "$url" | gunzip > "$dest"
    fi
}

# --- Run Prodigal on a genome ---
run_prodigal() {
    local fasta="$1" output="$2"
    prodigal -i "$fasta" -f gff -p single -o "$output" 2>/dev/null
}

# --- Header ---
echo "============================================================"
echo "  Prokrustes Benchmark Suite"
echo "============================================================"
if $RUN_PRODIGAL; then
    echo "  Baseline: $PRODIGAL_BIN"
fi
echo ""

# --- Summary files ---
SUMMARY="$RESULTS_DIR/benchmark_summary.tsv"
echo -e "Genome\tSpecies\tSize_Mb\tTool\tRef_CDS\tPred_CDS\tTP\tFP\tFN\tPrecision\tRecall\tF1\tStart_Acc" > "$SUMMARY"

for entry in "${GENOMES[@]}"; do
    IFS='|' read -r key label fasta_url gff_url <<< "$entry"

    if [ -n "$FILTER" ] && [ "$FILTER" != "$key" ]; then
        continue
    fi

    echo "--- $label ---"

    # Download
    download "$fasta_url" "$DATA_DIR/${key}.fasta"
    download "$gff_url" "$DATA_DIR/${key}.gff"

    GENOME_BYTES=$(wc -c < "$DATA_DIR/${key}.fasta" | tr -d ' ')
    SIZE_MB=$(echo "scale=1; $GENOME_BYTES / 1048576" | bc)

    # --- Prokrustes ---
    echo "  [Prokrustes] Annotating ($SIZE_MB Mb)..."
    START_TIME=$SECONDS
    "$BINARY" "$DATA_DIR/${key}.fasta" > "$RESULTS_DIR/${key}_prokrustes.gff" 2> "$RESULTS_DIR/${key}_prokrustes_stderr.log"
    ELAPSED=$((SECONDS - START_TIME))

    EVAL_PK=$(python3 "$SCRIPT_DIR/evaluate.py" "$RESULTS_DIR/${key}_prokrustes.gff" "$DATA_DIR/${key}.gff" --tsv)
    F1_PK=$(echo "$EVAL_PK" | cut -f8)
    echo "  [Prokrustes] F1=$F1_PK  (${ELAPSED}s)"
    echo -e "${key}\t${label}\t${SIZE_MB}\tProkrustes\t${EVAL_PK}" >> "$SUMMARY"

    # --- Prodigal ---
    if $RUN_PRODIGAL; then
        echo "  [Prodigal]   Annotating..."
        START_TIME=$SECONDS
        run_prodigal "$DATA_DIR/${key}.fasta" "$RESULTS_DIR/${key}_prodigal.gff"
        ELAPSED=$((SECONDS - START_TIME))

        EVAL_PD=$(python3 "$SCRIPT_DIR/evaluate.py" "$RESULTS_DIR/${key}_prodigal.gff" "$DATA_DIR/${key}.gff" --tsv)
        F1_PD=$(echo "$EVAL_PD" | cut -f8)
        echo "  [Prodigal]   F1=$F1_PD  (${ELAPSED}s)"
        echo -e "${key}\t${label}\t${SIZE_MB}\tProdigal\t${EVAL_PD}" >> "$SUMMARY"
    fi

    echo ""
done

# --- Print summary table ---
echo "============================================================"
echo "  Summary"
echo "============================================================"
echo ""

# Pretty-print comparison
python3 << 'PYEOF'
import sys

results = {}
with open("RESULTS_DIR/benchmark_summary.tsv".replace("RESULTS_DIR", sys.argv[1])) as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 13:
            continue
        key, species, size, tool = parts[0], parts[1], parts[2], parts[3]
        f1, prec, recall = parts[11], parts[9], parts[10]
        start_acc = parts[12]
        if key not in results:
            results[key] = {"species": species, "size": size}
        results[key][tool] = {"f1": f1, "prec": prec, "recall": recall, "start_acc": start_acc}

print(f"{'Genome':<28} {'Size':>5}  {'Prokrustes F1':>14}  {'Prodigal F1':>12}  {'Delta':>7}")
print("-" * 75)

for key, data in results.items():
    pk = data.get("Prokrustes", {})
    pd = data.get("Prodigal", {})
    pk_f1 = pk.get("f1", "-")
    pd_f1 = pd.get("f1", "-")

    delta = ""
    if pk_f1 != "-" and pd_f1 != "-":
        d = float(pk_f1) - float(pd_f1)
        delta = f"{d:+.4f}"

    print(f"{data['species'][:28]:<28} {data['size']:>5}  {pk_f1:>14}  {pd_f1:>12}  {delta:>7}")

print()
PYEOF
"$RESULTS_DIR"

echo "Full results: $RESULTS_DIR/"
echo "Summary TSV:  $SUMMARY"
