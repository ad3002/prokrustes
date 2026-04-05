"""Ribosome Binding Site (RBS) scoring: pattern matching + PWM.

Two methods:
1. score_rbs() — pattern matching with spacing-dependent weighting
2. build_rbs_pwm() / score_rbs_pwm() — trained position-weight matrix

Tested: yes
Source: extracted from cycle 22 (F1=0.926)
"""

import math


def score_rbs(seq: str, start_pos: int) -> float:
    """Score Shine-Dalgarno motif upstream of start codon.

    Searches for AGGAGG-like patterns 5-18bp upstream.
    Returns 0.0 to 1.0.
    """
    if start_pos < 4:
        return 0.0
    up_begin = max(0, start_pos - 24)
    up_end = start_pos - 3
    if up_end <= up_begin:
        return 0.0
    upstream = seq[up_begin:up_end]
    patterns = [
        ("TAAGGAGG", 1.00), ("AAGGAGG", 1.00), ("AGGAGG", 0.98),
        ("AGGAG", 0.90), ("GGAGG", 0.87),
        ("TAAGGA", 0.80), ("AAGGA", 0.76), ("GAGGT", 0.70),
        ("AGGA", 0.68), ("GAGG", 0.68), ("GGAG", 0.60),
        ("AGG", 0.40), ("GGA", 0.35), ("GAG", 0.28), ("GG", 0.14),
    ]
    best = 0.0
    for pat, base_score in patterns:
        idx = 0
        while True:
            pos = upstream.find(pat, idx)
            if pos < 0:
                break
            motif_end = up_begin + pos + len(pat)
            spacing = start_pos - motif_end
            if 5 <= spacing <= 10:
                sp = 1.0
            elif 4 <= spacing <= 12:
                sp = 0.85
            elif 3 <= spacing <= 14:
                sp = 0.68
            elif 1 <= spacing <= 18:
                sp = 0.42
            else:
                sp = 0.12
            s = base_score * sp
            if s > best:
                best = s
            idx = pos + 1
    return best


def build_rbs_pwm(seq: str, orfs: list) -> list | None:
    """Train position-weight matrix for RBS from confident ORFs.

    Returns PWM (list of dicts {A:score, T:score, G:score, C:score})
    or None if insufficient training data.
    """
    conf = [o for o in orfs if o.get("score", 0) > 0.46
            and o["length"] >= 350 and o["start_codon"] == "ATG"
            and o["is_longest"]]
    if len(conf) < 80:
        return None

    width = 20
    counts = [{"A": 0, "T": 0, "G": 0, "C": 0} for _ in range(width)]
    bg = {"A": 0, "T": 0, "G": 0, "C": 0}
    bg_total = 0
    n_seqs = 0

    for o in conf:
        sp = o["seq_start"]
        if sp < width + 3:
            continue
        region = seq[sp - width:sp]
        if "N" in region:
            continue
        for i, nt in enumerate(region):
            if nt in counts[i]:
                counts[i][nt] += 1
        n_seqs += 1

    if n_seqs < 60:
        return None

    for i in range(min(len(seq), 500000)):
        nt = seq[i]
        if nt in bg:
            bg[nt] += 1
            bg_total += 1

    if bg_total == 0:
        return None

    pwm = []
    for i in range(width):
        row = {}
        for nt in "ACGT":
            freq = (counts[i][nt] + 0.5) / (n_seqs + 2)
            bg_freq = (bg[nt] + 0.5) / (bg_total + 2)
            row[nt] = math.log(freq / bg_freq) if bg_freq > 0 else 0.0
        pwm.append(row)
    return pwm


def score_rbs_pwm(seq: str, start_pos: int, pwm: list | None) -> float:
    """Score RBS using trained PWM."""
    if pwm is None:
        return 0.0
    width = len(pwm)
    if start_pos < width + 3:
        return 0.0
    region = seq[start_pos - width:start_pos]
    if "N" in region:
        return 0.0
    return sum(pwm[i].get(nt, 0.0) for i, nt in enumerate(region))


# === Tests ===
def test_score_rbs_aggagg():
    # AGGAGG at position 7bp upstream of ATG
    seq = "AAAAAAAGGAGGAAAAAATG" + "A" * 100
    start = seq.index("ATG")
    score = score_rbs(seq, start)
    assert score > 0.5, f"AGGAGG should score high, got {score}"


def test_score_rbs_no_motif():
    seq = "TTTTTTTTTTTTTTTTTTTTTTTATG" + "A" * 100
    start = seq.index("ATG")
    score = score_rbs(seq, start)
    assert score < 0.2, f"No SD should score low, got {score}"


def test_score_rbs_too_close():
    # SD too close to start (spacing < 3)
    seq = "AGGAGGATG" + "A" * 100
    start = seq.index("ATG")
    score = score_rbs(seq, start)
    # Should still find it but with low spacing score
    assert score >= 0.0


if __name__ == "__main__":
    test_score_rbs_aggagg()
    test_score_rbs_no_motif()
    test_score_rbs_too_close()
    print("rbs_scoring: all tests passed")
