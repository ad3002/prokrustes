"""ORF finder: scan all 3 reading frames for open reading frames.

Finds all ORFs bounded by start codons (ATG/GTG/TTG) and stop codons
(TAA/TAG/TGA) in a single strand. Returns list of ORF dicts with
genomic coordinates.

Tested: yes
Source: extracted from cycle 22 (F1=0.926)
"""

START_CODONS = frozenset({"ATG", "GTG", "TTG"})
STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
MIN_GENE_LEN = 90


def find_all_orfs(seq: str, strand: str, genome_len: int, min_len: int = MIN_GENE_LEN) -> list:
    """Find all ORFs in seq (one strand).

    Args:
        seq: DNA sequence (forward or reverse complement)
        strand: '+' or '-'
        genome_len: total genome length (for coordinate conversion)
        min_len: minimum ORF length in bp

    Returns:
        List of ORF dicts with keys:
        start, end, strand, length, seq_start, seq_end,
        is_longest, start_codon, frame, stop_pos
    """
    orfs = []
    n = len(seq)
    for frame in range(3):
        # Collect stop codon positions in this frame
        stops = []
        for i in range(frame, n - 2, 3):
            if seq[i:i+3] in STOP_CODONS:
                stops.append(i)

        # Between consecutive stops, find all start codons
        prev_stop_end = frame
        for stop_pos in stops:
            stop_end = stop_pos + 3
            starts = []
            for j in range(prev_stop_end, stop_pos, 3):
                if seq[j:j+3] in START_CODONS:
                    starts.append(j)

            for idx, start_pos in enumerate(starts):
                orf_len = stop_end - start_pos
                if orf_len < min_len:
                    continue
                if strand == "+":
                    g_start = start_pos + 1
                    g_end = stop_end
                else:
                    g_start = genome_len - stop_end + 1
                    g_end = genome_len - start_pos

                orfs.append({
                    "start": g_start, "end": g_end, "strand": strand,
                    "length": orf_len, "seq_start": start_pos, "seq_end": stop_end,
                    "is_longest": idx == 0,
                    "start_codon": seq[start_pos:start_pos+3],
                    "frame": frame, "stop_pos": stop_pos,
                })
            prev_stop_end = stop_end
    return orfs


# === Tests ===
def test_simple_orf():
    # ATG...TAA in frame 0
    seq = "ATGAAAAAATAA"  # ATG + AAA + AAA + TAA = 12bp
    orfs = find_all_orfs(seq, "+", len(seq), min_len=9)
    assert len(orfs) >= 1
    assert orfs[0]["start_codon"] == "ATG"
    assert orfs[0]["length"] == 12


def test_no_orf_too_short():
    seq = "ATGTAA"  # 6bp < 90bp minimum
    orfs = find_all_orfs(seq, "+", len(seq))
    assert len(orfs) == 0


def test_multiple_frames():
    # ORFs can be in different frames
    seq = "A" + "ATG" + "AAA" * 40 + "TAA"  # frame 1, 126bp ORF
    orfs = find_all_orfs(seq, "+", len(seq))
    assert any(o["frame"] == 1 for o in orfs)


def test_minus_strand_coordinates():
    seq = "ATG" + "AAA" * 40 + "TAA"  # 126bp
    genome_len = 1000
    orfs = find_all_orfs(seq, "-", genome_len)
    assert len(orfs) >= 1
    # Minus strand: coordinates are flipped
    assert orfs[0]["strand"] == "-"


if __name__ == "__main__":
    test_simple_orf()
    test_no_orf_too_short()
    test_multiple_frames()
    test_minus_strand_coordinates()
    print("orf_finder: all tests passed")
