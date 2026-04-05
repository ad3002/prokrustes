#!/usr/bin/env python3
"""Extract rRNA and tRNA regions from NCBI GFF to use as ncRNA mask."""
import sys

def extract_ncrna(gff_path, output_path):
    with open(gff_path) as f, open(output_path, 'w') as out:
        out.write("##gff-version 3\n")
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            if fields[2] in ('rRNA', 'tRNA'):
                out.write(line)

if __name__ == "__main__":
    extract_ncrna(sys.argv[1], sys.argv[2])
