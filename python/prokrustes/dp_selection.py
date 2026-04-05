"""Dynamic programming gene selection: optimal non-overlapping gene tiling.

Instead of greedy overlap resolution, finds globally optimal set of genes
that maximizes total score while respecting overlap constraints.

Tested: yes
Source: extracted from cycle 22 (F1=0.926)
"""


def dp_gene_selection(orfs: list, max_overlap_bp: int = 60) -> list:
    """Select optimal non-overlapping genes via weighted interval scheduling.

    Each ORF has a score. Find subset with max total score such that
    no two ORFs overlap more than max_overlap_bp.

    Args:
        orfs: list of ORF dicts with 'start', 'end', 'score' keys
        max_overlap_bp: maximum allowed overlap between selected genes

    Returns:
        Selected ORFs.
    """
    if not orfs:
        return []

    # Sort by end position
    sorted_orfs = sorted(orfs, key=lambda o: o["end"])

    n = len(sorted_orfs)
    dp = [0.0] * (n + 1)
    parent = [-1] * (n + 1)

    for i in range(1, n + 1):
        orf = sorted_orfs[i - 1]
        w = orf.get("score", orf.get("total_score", 0.0))

        # Find latest compatible ORF
        j = _find_latest_compatible(sorted_orfs, i - 1, max_overlap_bp)

        include = dp[j + 1] + w
        exclude = dp[i - 1]

        if include > exclude:
            dp[i] = include
            parent[i] = j + 1
        else:
            dp[i] = exclude
            parent[i] = -(i - 1)  # negative = exclude

    # Backtrack
    selected = []
    i = n
    while i > 0:
        orf = sorted_orfs[i - 1]
        w = orf.get("score", orf.get("total_score", 0.0))
        j = _find_latest_compatible(sorted_orfs, i - 1, max_overlap_bp)
        if dp[j + 1] + w > dp[i - 1]:
            selected.append(orf)
            i = j + 1
        else:
            i -= 1

    return selected


def _find_latest_compatible(orfs: list, idx: int, max_overlap: int) -> int:
    """Binary search for latest ORF that doesn't overlap too much with orfs[idx]."""
    target_start = orfs[idx]["start"] + max_overlap
    lo, hi = 0, idx - 1
    result = -1
    while lo <= hi:
        mid = (lo + hi) // 2
        if orfs[mid]["end"] <= target_start:
            result = mid
            lo = mid + 1
        else:
            hi = mid - 1
    return result


def iterative_connection_dp(candidates: list, max_overlap_bp: int = 60,
                            n_iterations: int = 4) -> list:
    """Multi-pass DP: each iteration boosts adjacency scores.

    After first DP pass, genes that are adjacent (close together)
    get bonus scores. Re-run DP to pick up operons.
    """
    current = candidates[:]
    selected = dp_gene_selection(current, max_overlap_bp)

    for iteration in range(1, n_iterations):
        # Boost genes adjacent to selected ones
        selected_set = {(g["start"], g["end"]) for g in selected}
        selected_ends = sorted(g["end"] for g in selected)

        for orf in current:
            if (orf["start"], orf["end"]) in selected_set:
                continue
            # Check if adjacent to a selected gene
            for se in selected_ends:
                gap = orf["start"] - se
                if 0 <= gap <= 300:
                    orf["score"] = orf.get("score", 0) * 1.15
                    break
                if gap > 300:
                    break

        selected = dp_gene_selection(current, max_overlap_bp)

    return selected


# === Tests ===
def test_dp_no_overlap():
    orfs = [
        {"start": 1, "end": 100, "score": 10},
        {"start": 200, "end": 300, "score": 20},
        {"start": 400, "end": 500, "score": 15},
    ]
    selected = dp_gene_selection(orfs)
    assert len(selected) == 3  # no overlaps, all selected


def test_dp_with_overlap():
    orfs = [
        {"start": 1, "end": 200, "score": 10},
        {"start": 150, "end": 300, "score": 20},  # overlaps first
    ]
    selected = dp_gene_selection(orfs, max_overlap_bp=0)
    assert len(selected) == 1
    assert selected[0]["score"] == 20  # higher score wins


def test_dp_allows_small_overlap():
    orfs = [
        {"start": 1, "end": 200, "score": 10},
        {"start": 198, "end": 400, "score": 10},  # 2bp overlap
    ]
    selected = dp_gene_selection(orfs, max_overlap_bp=4)
    assert len(selected) == 2  # small overlap allowed


if __name__ == "__main__":
    test_dp_no_overlap()
    test_dp_with_overlap()
    test_dp_allows_small_overlap()
    print("dp_selection: all tests passed")
