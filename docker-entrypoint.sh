#!/usr/bin/env bash
set -euo pipefail

case "${1:-}" in
    benchmark)
        shift
        exec bash /opt/prokrustes/scripts/benchmark.sh "$@"
        ;;
    annotate)
        shift
        exec prokrustes "$@"
        ;;
    --help|-h|"")
        echo "Prokrustes — prokaryotic genome annotation"
        echo ""
        echo "Usage:"
        echo "  docker run prokrustes benchmark                    # full 10-genome benchmark vs Prodigal"
        echo "  docker run prokrustes benchmark ecoli_k12          # single genome"
        echo "  docker run prokrustes benchmark --no-prodigal      # without Prodigal baseline"
        echo "  docker run prokrustes annotate /data/ecoli_k12.fasta   # annotate a genome"
        echo "  docker run prokrustes annotate your_genome.fasta       # annotate custom genome"
        echo ""
        echo "Pre-loaded genomes in /data/:"
        ls /data/*.fasta 2>/dev/null | while read f; do
            name=$(basename "$f" .fasta)
            size=$(wc -c < "$f" | tr -d ' ')
            mb=$(echo "scale=1; $size / 1048576" | bc)
            echo "  $name (${mb} Mb)"
        done
        echo ""
        echo "Mount your own genomes:"
        echo "  docker run -v /path/to/genomes:/workspace prokrustes annotate /workspace/genome.fasta"
        ;;
    *)
        # Pass through to prokrustes binary
        exec prokrustes "$@"
        ;;
esac
