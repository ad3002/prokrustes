"""I/O utilities: read FASTA, reverse complement.

Tested: yes
Source: extracted from cycle 22 (F1=0.926)
"""

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_fasta(path: str) -> str:
    """Read FASTA file, return uppercase sequence string."""
    parts = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                parts.append(line.strip().upper())
    return "".join(parts)


def reverse_complement(seq: str) -> str:
    """Reverse complement of DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


# === Tests ===
def test_read_fasta(tmp_path=None):
    import tempfile, os
    fd, path = tempfile.mkstemp(suffix=".fasta")
    with os.fdopen(fd, 'w') as f:
        f.write(">test\nATGC\nGGCC\n")
    assert read_fasta(path) == "ATGCGGCC"
    os.unlink(path)


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAAAA") == "TTTTTT"
    assert reverse_complement("") == ""


if __name__ == "__main__":
    test_read_fasta()
    test_reverse_complement()
    print("io_utils: all tests passed")
