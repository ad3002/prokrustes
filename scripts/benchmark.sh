#!/usr/bin/env bash
set -euo pipefail

# Prokrustes — parallel benchmark with Prodigal baseline
# Runs all genomes in parallel (one per core), then prints comparison table.
#
# Usage:
#   ./scripts/benchmark.sh                  # all genomes, parallel
#   ./scripts/benchmark.sh ecoli_k12        # single genome
#   ./scripts/benchmark.sh --no-prodigal    # skip Prodigal baseline

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# Auto-detect Docker vs local environment
if [ -f /usr/local/bin/prokrustes ] && [ -d /data ]; then
    DATA_DIR="/data"
    RESULTS_DIR="/results"
    BINARY="/usr/local/bin/prokrustes"
else
    DATA_DIR="$PROJECT_DIR/data"
    RESULTS_DIR="$PROJECT_DIR/results"
    BINARY="$PROJECT_DIR/target/release/prokrustes"
fi

mkdir -p "$DATA_DIR" "$RESULTS_DIR"

# --- Genome registry ---
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
if $RUN_PRODIGAL && ! command -v prodigal &>/dev/null; then
    echo "NOTE: Prodigal not found. Install: conda install -c bioconda prodigal"
    echo "      Running without baseline comparison."
    echo ""
    RUN_PRODIGAL=false
fi

# --- Build if needed ---
if [ ! -f "$BINARY" ]; then
    echo "Building prokrustes (release)..."
    cd "$PROJECT_DIR" && cargo build --release 2>&1 | tail -3
fi

download() {
    local url="$1" dest="$2"
    if [ ! -f "$dest" ]; then
        curl -sL "$url" | gunzip > "$dest"
    fi
}

# --- Download all data first (sequential, fast) ---
echo "Downloading reference genomes..."
for entry in "${GENOMES[@]}"; do
    IFS='|' read -r key label fasta_url gff_url <<< "$entry"
    if [ -n "$FILTER" ] && [ "$FILTER" != "$key" ]; then continue; fi
    download "$fasta_url" "$DATA_DIR/${key}.fasta"
    download "$gff_url" "$DATA_DIR/${key}.gff"
done
echo "Done."
echo ""

# --- Worker function: annotate + evaluate one genome ---
run_one() {
    local key="$1" label="$2"
    local log="$RESULTS_DIR/${key}.log"

    local genome_bytes
    genome_bytes=$(wc -c < "$DATA_DIR/${key}.fasta" | tr -d ' ')
    local size_mb
    size_mb=$(echo "scale=1; $genome_bytes / 1048576" | bc)

    # ncRNA mask
    local ncrna_flag=""
    if command -v barrnap &>/dev/null; then
        local ncrna_gff="$DATA_DIR/${key}_ncrna.gff"
        if [ ! -f "$ncrna_gff" ]; then
            barrnap "$DATA_DIR/${key}.fasta" > "$ncrna_gff" 2>/dev/null || true
        fi
        if [ -s "$ncrna_gff" ]; then
            ncrna_flag="--ncrna $ncrna_gff"
        fi
    fi

    # Prokrustes
    local pk_start=$SECONDS
    $BINARY "$DATA_DIR/${key}.fasta" $ncrna_flag > "$RESULTS_DIR/${key}_prokrustes.gff" 2> "$RESULTS_DIR/${key}_prokrustes_stderr.log"
    local pk_time=$((SECONDS - pk_start))

    local eval_pk
    eval_pk=$(python3 "$SCRIPT_DIR/evaluate.py" "$RESULTS_DIR/${key}_prokrustes.gff" "$DATA_DIR/${key}.gff" --tsv)
    local f1_pk
    f1_pk=$(echo "$eval_pk" | cut -f8)

    # Write per-genome result
    echo -e "${key}\t${label}\t${size_mb}\tProkrustes\t${eval_pk}\t${pk_time}" >> "$RESULTS_DIR/${key}_results.tsv"

    # Prodigal
    if $RUN_PRODIGAL; then
        local pd_start=$SECONDS
        prodigal -i "$DATA_DIR/${key}.fasta" -f gff -p single -o "$RESULTS_DIR/${key}_prodigal.gff" 2>/dev/null
        local pd_time=$((SECONDS - pd_start))

        local eval_pd
        eval_pd=$(python3 "$SCRIPT_DIR/evaluate.py" "$RESULTS_DIR/${key}_prodigal.gff" "$DATA_DIR/${key}.gff" --tsv)
        local f1_pd
        f1_pd=$(echo "$eval_pd" | cut -f8)

        echo -e "${key}\t${label}\t${size_mb}\tProdigal\t${eval_pd}\t${pd_time}" >> "$RESULTS_DIR/${key}_results.tsv"
        echo "  done: $label — Prokrustes F1=$f1_pk (${pk_time}s) | Prodigal F1=$f1_pd (${pd_time}s)"
    else
        echo "  done: $label — Prokrustes F1=$f1_pk (${pk_time}s)"
    fi
}

# --- Header ---
NCORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
echo "============================================================"
echo "  Prokrustes Benchmark Suite (parallel, $NCORES cores)"
echo "============================================================"
if $RUN_PRODIGAL; then
    echo "  Baseline: prodigal -p single"
fi
echo ""

# Clean per-genome result files
for entry in "${GENOMES[@]}"; do
    IFS='|' read -r key label fasta_url gff_url <<< "$entry"
    rm -f "$RESULTS_DIR/${key}_results.tsv"
done

# --- Launch all genomes in parallel ---
TOTAL_START=$SECONDS
PIDS=()

for entry in "${GENOMES[@]}"; do
    IFS='|' read -r key label fasta_url gff_url <<< "$entry"
    if [ -n "$FILTER" ] && [ "$FILTER" != "$key" ]; then continue; fi

    run_one "$key" "$label" &
    PIDS+=($!)

    # Limit parallelism to NCORES
    while [ ${#PIDS[@]} -ge "$NCORES" ]; do
        # Wait for any one to finish
        for i in "${!PIDS[@]}"; do
            if ! kill -0 "${PIDS[$i]}" 2>/dev/null; then
                wait "${PIDS[$i]}" || true
                unset 'PIDS[$i]'
            fi
        done
        PIDS=("${PIDS[@]}")
        sleep 0.1
    done
done

# Wait for remaining
for pid in "${PIDS[@]}"; do
    wait "$pid" || true
done

TOTAL_ELAPSED=$((SECONDS - TOTAL_START))

# --- Merge results ---
SUMMARY="$RESULTS_DIR/benchmark_summary.tsv"
echo -e "Genome\tSpecies\tSize_Mb\tTool\tRef_CDS\tPred_CDS\tTP\tFP\tFN\tPrecision\tRecall\tF1\tStart_Acc\tTime_s" > "$SUMMARY"

for entry in "${GENOMES[@]}"; do
    IFS='|' read -r key label fasta_url gff_url <<< "$entry"
    if [ -n "$FILTER" ] && [ "$FILTER" != "$key" ]; then continue; fi
    if [ -f "$RESULTS_DIR/${key}_results.tsv" ]; then
        cat "$RESULTS_DIR/${key}_results.tsv" >> "$SUMMARY"
    fi
done

# --- Print summary table ---
echo ""
echo "============================================================"
echo "  Results (total wall time: ${TOTAL_ELAPSED}s)"
echo "============================================================"
echo ""

BENCH_RESULTS="$RESULTS_DIR" python3 << 'PYEOF'
import os

results_dir = os.environ["BENCH_RESULTS"]
results = {}
order = []
with open(results_dir + "/benchmark_summary.tsv") as f:
    header = f.readline()
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) < 4:
            continue
        key, species, size, tool = parts[0], parts[1], parts[2], parts[3]
        f1 = parts[11] if len(parts) > 11 else "-"
        time_s = parts[13] if len(parts) > 13 else "-"
        if key not in results:
            results[key] = {"species": species, "size": size}
            order.append(key)
        results[key][tool] = {"f1": f1, "time": time_s}

print(f"{'Genome':<30} {'Size':>5}  {'Prokrustes':>11} {'t':>4}  {'Prodigal':>9} {'t':>4}  {'Delta':>7}")
print("-" * 80)

for key in order:
    data = results[key]
    pk = data.get("Prokrustes", {})
    pd = data.get("Prodigal", {})
    pk_f1 = pk.get("f1", "-")
    pd_f1 = pd.get("f1", "-")
    pk_t = pk.get("time", "-")
    pd_t = pd.get("time", "-")

    delta = ""
    try:
        if pk_f1 != "-" and pd_f1 != "-":
            delta = f"{float(pk_f1) - float(pd_f1):+.4f}"
    except ValueError:
        pass

    pk_str = f"F1={pk_f1}" if pk_f1 != "-" else "-"
    pd_str = f"F1={pd_f1}" if pd_f1 != "-" else "-"
    pk_t_str = f"{pk_t}s" if pk_t != "-" else "-"
    pd_t_str = f"{pd_t}s" if pd_t != "-" else "-"

    print(f"{data['species'][:30]:<30} {data['size']:>5}  {pk_str:>11} {pk_t_str:>4}  {pd_str:>9} {pd_t_str:>4}  {delta:>7}")

print()
PYEOF

echo "Full results: $RESULTS_DIR/"
echo "Summary TSV:  $SUMMARY"
