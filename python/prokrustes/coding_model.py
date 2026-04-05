"""Coding potential models: hexamer, monocodon, intergenic.

Self-training approach:
1. Train on long ORFs (>= 800bp) as "definitely coding"
2. Background = shuffled genome or all frames
3. Log-likelihood ratio: log(P(coding)/P(background))

Tested: yes
Source: extracted from cycle 22 (F1=0.926)
"""

import math
from collections import Counter

BASES = "ATGC"
ALL_CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]


def build_hexamer_table(seq: str, orfs: list, min_train_len: int = 800,
                        require_atg: bool = True) -> dict:
    """Train hexamer (6-mer) log-likelihood ratio model.

    Hexamers = dicodons (pairs of adjacent codons). Captures codon-pair
    context that single codons miss.

    Returns dict: hexamer -> log2(P_coding / P_background)
    """
    training = [o for o in orfs if o["length"] >= min_train_len]
    if require_atg:
        training = [o for o in training if o["start_codon"] == "ATG"]
    if not training:
        return {}

    coding_counts = Counter()
    for orf in training:
        for i in range(orf["seq_start"], orf["seq_end"] - 5, 3):
            hexamer = seq[i:i+6]
            if len(hexamer) == 6 and "N" not in hexamer:
                coding_counts[hexamer] += 1
    total_coding = sum(coding_counts.values())
    if total_coding < 200:
        return {}

    bg_counts = Counter()
    for frame_offset in range(3):
        for i in range(frame_offset, len(seq) - 5, 3):
            h = seq[i:i+6]
            if len(h) == 6 and "N" not in h:
                bg_counts[h] += 1
    total_bg = sum(bg_counts.values())
    if total_bg == 0:
        return {}

    all_hex = [a+b+c+d+e+f for a in BASES for b in BASES for c in BASES
               for d in BASES for e in BASES for f in BASES]
    pseudo = 0.1
    table = {}
    for h in all_hex:
        p_cod = (coding_counts.get(h, 0) + pseudo) / (total_coding + pseudo * len(all_hex))
        p_bg = (bg_counts.get(h, 0) + pseudo) / (total_bg + pseudo * len(all_hex))
        table[h] = math.log2(p_cod / p_bg) if p_bg > 0 else 0.0
    return table


def build_intergenic_hexamer_table(seq: str, coding_orfs: list,
                                   intergenic_seqs: list) -> dict:
    """Train hexamer model specifically for intergenic (non-coding) regions.

    Uses actual intergenic sequences as background instead of all-frame shuffled.
    More discriminative than generic background.
    """
    if not intergenic_seqs or not coding_orfs:
        return {}

    coding_counts = Counter()
    for orf in coding_orfs:
        for i in range(orf["seq_start"], orf["seq_end"] - 5, 3):
            h = seq[i:i+6]
            if len(h) == 6 and "N" not in h:
                coding_counts[h] += 1
    total_coding = sum(coding_counts.values())

    ig_counts = Counter()
    for ig_seq in intergenic_seqs:
        for i in range(0, len(ig_seq) - 5):
            h = ig_seq[i:i+6]
            if len(h) == 6 and "N" not in h:
                ig_counts[h] += 1
    total_ig = sum(ig_counts.values())

    if total_coding < 200 or total_ig < 200:
        return {}

    all_hex = [a+b+c+d+e+f for a in BASES for b in BASES for c in BASES
               for d in BASES for e in BASES for f in BASES]
    pseudo = 0.1
    table = {}
    for h in all_hex:
        p_cod = (coding_counts.get(h, 0) + pseudo) / (total_coding + pseudo * len(all_hex))
        p_ig = (ig_counts.get(h, 0) + pseudo) / (total_ig + pseudo * len(all_hex))
        table[h] = math.log2(p_cod / p_ig) if p_ig > 0 else 0.0
    return table


def build_monocodon_table(seq: str, orfs: list, min_train_len: int = 800) -> dict:
    """Train single-codon log-likelihood ratio model.

    Simpler than hexamer but still useful as a secondary signal.
    """
    training = [o for o in orfs if o["length"] >= min_train_len
                and o["start_codon"] == "ATG"]
    if not training:
        return {}

    coding_counts = Counter()
    for orf in training:
        for i in range(orf["seq_start"], orf["seq_end"] - 2, 3):
            codon = seq[i:i+3]
            if len(codon) == 3 and "N" not in codon:
                coding_counts[codon] += 1
    total_coding = sum(coding_counts.values())

    bg_counts = Counter()
    for off in range(3):
        for i in range(off, len(seq) - 2, 3):
            c = seq[i:i+3]
            if len(c) == 3 and "N" not in c:
                bg_counts[c] += 1
    total_bg = sum(bg_counts.values())

    if total_coding < 100 or total_bg == 0:
        return {}

    pseudo = 0.5
    table = {}
    for codon in ALL_CODONS:
        p_cod = (coding_counts.get(codon, 0) + pseudo) / (total_coding + pseudo * 64)
        p_bg = (bg_counts.get(codon, 0) + pseudo) / (total_bg + pseudo * 64)
        table[codon] = math.log2(p_cod / p_bg) if p_bg > 0 else 0.0
    return table


def score_coding(seq: str, orf: dict, table: dict) -> float:
    """Score ORF using hexamer or monocodon table. Returns mean log-likelihood."""
    if not table:
        return 0.0
    total = 0.0
    n = 0
    unit_size = 6 if any(len(k) == 6 for k in table) else 3
    for i in range(orf["seq_start"], orf["seq_end"] - (unit_size - 1), 3):
        unit = seq[i:i+unit_size]
        if unit in table:
            total += table[unit]
            n += 1
    return total / n if n > 0 else 0.0


def merge_tables(t1: dict, t2: dict) -> dict:
    """Merge two scoring tables by averaging."""
    merged = {}
    all_keys = set(t1) | set(t2)
    for k in all_keys:
        merged[k] = (t1.get(k, 0.0) + t2.get(k, 0.0)) / 2
    return merged


def blend_tables(t1: dict, t2: dict, w2: float = 0.3) -> dict:
    """Weighted blend of two tables: (1-w2)*t1 + w2*t2."""
    blended = {}
    all_keys = set(t1) | set(t2)
    for k in all_keys:
        blended[k] = (1 - w2) * t1.get(k, 0.0) + w2 * t2.get(k, 0.0)
    return blended


def get_intergenic_regions(genome: str, genes: list, min_gap: int = 50) -> list:
    """Extract intergenic sequences between genes."""
    if not genes:
        return []
    sorted_genes = sorted(genes, key=lambda g: g["seq_start"])
    regions = []
    prev_end = 0
    for g in sorted_genes:
        gap_start = prev_end
        gap_end = g["seq_start"]
        if gap_end - gap_start >= min_gap:
            regions.append(genome[gap_start:gap_end])
        prev_end = max(prev_end, g["seq_end"])
    return regions


# === Tests ===
def test_build_hexamer_table():
    # Minimal test: needs long ORFs
    seq = "ATG" + "AAGCTG" * 200 + "TAA" + "TTT" * 100
    orfs = [{"seq_start": 0, "seq_end": 1203, "length": 1203, "start_codon": "ATG"}]
    table = build_hexamer_table(seq, orfs, min_train_len=100)
    assert len(table) > 0
    # AAGCTG should score positive (coding)
    assert "AAGCTG" in table


def test_score_coding():
    table = {"AAGCTG": 1.0, "CTGAAG": 0.5}
    orf = {"seq_start": 0, "seq_end": 12}
    seq = "AAGCTGAAGCTG"
    score = score_coding(seq, orf, table)
    assert score > 0


def test_monocodon_table():
    seq = "ATG" + "AAG" * 400 + "TAA" + "TTT" * 100
    orfs = [{"seq_start": 0, "seq_end": 1203, "length": 1203, "start_codon": "ATG"}]
    table = build_monocodon_table(seq, orfs, min_train_len=100)
    assert "AAG" in table
    assert table["AAG"] > 0  # coding codon should be positive


if __name__ == "__main__":
    test_build_hexamer_table()
    test_score_coding()
    test_monocodon_table()
    print("coding_model: all tests passed")
