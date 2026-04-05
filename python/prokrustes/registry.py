"""Registry of available library parts.

Each part has:
- name: unique identifier
- module: Python module path
- description: what it does
- status: tested/untested/broken
- f1_contribution: how much it helps (estimated)
"""

PARTS = {
    "io_utils": {
        "module": "gene_finder_lib.io_utils",
        "functions": ["read_fasta", "reverse_complement"],
        "description": "FASTA I/O and reverse complement",
        "status": "tested",
    },
    "orf_finder": {
        "module": "gene_finder_lib.orf_finder",
        "functions": ["find_all_orfs"],
        "description": "Find ORFs in all 3 reading frames of one strand",
        "status": "tested",
    },
    "rbs_scoring": {
        "module": "gene_finder_lib.rbs_scoring",
        "functions": ["score_rbs", "build_rbs_pwm", "score_rbs_pwm"],
        "description": "Shine-Dalgarno pattern matching + trained PWM",
        "status": "tested",
    },
    "coding_model": {
        "module": "gene_finder_lib.coding_model",
        "functions": [
            "build_hexamer_table", "build_intergenic_hexamer_table",
            "build_monocodon_table", "score_coding",
            "merge_tables", "blend_tables", "get_intergenic_regions",
        ],
        "description": "Self-training hexamer/monocodon log-likelihood models",
        "status": "tested",
    },
    "dp_selection": {
        "module": "gene_finder_lib.dp_selection",
        "functions": ["dp_gene_selection", "iterative_connection_dp"],
        "description": "Dynamic programming for optimal gene tiling",
        "status": "tested",
    },
}


def list_parts():
    """Print available parts."""
    for name, info in PARTS.items():
        print(f"  {name}: {info['description']} [{info['status']}]")
        for fn in info["functions"]:
            print(f"    - {fn}()")


def get_part_description(name: str) -> str:
    """Get description for LLM prompt."""
    if name not in PARTS:
        return f"Unknown part: {name}"
    p = PARTS[name]
    fns = ", ".join(p["functions"])
    return f"{name}: {p['description']}. Functions: {fns}. Status: {p['status']}."


def get_all_descriptions() -> str:
    """Get all part descriptions for LLM prompt."""
    lines = ["Available library parts:"]
    for name in PARTS:
        lines.append(f"  {get_part_description(name)}")
    return "\n".join(lines)


if __name__ == "__main__":
    list_parts()
