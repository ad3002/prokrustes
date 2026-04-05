"""Research pipeline that generates Rust gene finders.

LLM writes Rust code → cargo build → run binary → evaluate JSON output.
Evaluation stays in Python (fast, just JSON comparison).
"""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

WORK_DIR = Path("./annotator_work")
RUST_SRC = Path("./src/bin/gene_finder_candidate.rs")
BEST_RUST = Path("./src/bin/gene_finder.rs")


def call_llm(prompt: str, timeout: int = 300) -> str | None:
    try:
        result = subprocess.run(
            ["claude", "-p", prompt, "--output-format", "json",
             "--tools", "Read,Write", "--permission-mode", "bypassPermissions"],
            capture_output=True, text=True, timeout=timeout,
        )
        if result.returncode != 0:
            return None
        data = json.loads(result.stdout)
        return data.get("result", "")
    except Exception as e:
        print(f"  LLM error: {e}")
        return None


def extract_rust_code(response: str) -> str | None:
    if not response:
        return None
    for marker in ["```rust", "```json"]:
        if marker in response:
            after = response.split(marker, 1)[1]
            if "```" in after:
                block = after.split("```", 1)[0].strip()
                if marker == "```json":
                    try:
                        return json.loads(block).get("code")
                    except json.JSONDecodeError:
                        pass
                else:
                    return block
    try:
        return json.loads(response).get("code")
    except (json.JSONDecodeError, AttributeError):
        pass
    if "use std" in response:
        return response[response.index("use std"):]
    if "fn main" in response:
        return response[response.index("fn main"):]
    return None


def compile_rust(src_path: Path) -> tuple[bool, str]:
    """Compile Rust source, return (success, error_message)."""
    result = subprocess.run(
        ["cargo", "build", "--release", "--bin", "gene_finder_candidate"],
        capture_output=True, text=True, timeout=120,
    )
    if result.returncode == 0:
        return True, ""
    errors = [l for l in result.stderr.split('\n') if l.startswith("error")]
    return False, "\n".join(errors[:10])


def run_rust_binary(fasta_path: str, timeout: int = 120) -> tuple[str, float]:
    """Run compiled Rust binary, return (stdout, elapsed_seconds)."""
    start = time.time()
    result = subprocess.run(
        ["cargo", "run", "--release", "--bin", "gene_finder_candidate", "--", fasta_path],
        capture_output=True, text=True, timeout=timeout,
    )
    elapsed = time.time() - start
    if result.returncode != 0:
        raise RuntimeError(f"Runtime error: {result.stderr[:300]}")
    return result.stdout, elapsed


def evaluate_predictions(pred_json: str, gff_path: str) -> dict | None:
    try:
        predicted = json.loads(pred_json)
    except json.JSONDecodeError:
        if '[' in pred_json:
            try:
                predicted = json.loads(pred_json[pred_json.index('['):pred_json.rindex(']')+1])
            except (json.JSONDecodeError, ValueError):
                return None
        else:
            return None

    known = parse_gff_cds(gff_path)
    tp, fp, correct_starts = 0, 0, 0
    matched = [False] * len(known)

    for pred in predicted:
        best_ov, best_idx = 0.0, None
        for i, gene in enumerate(known):
            if matched[i] or pred.get("strand", "+") != gene["strand"]:
                continue
            ov_s = max(pred["start"], gene["start"])
            ov_e = min(pred["end"], gene["end"])
            if ov_s >= ov_e:
                continue
            ov_len = ov_e - ov_s
            pred_len = max(pred["end"] - pred["start"], 1)
            gene_len = max(gene["end"] - gene["start"], 1)
            recip = min(ov_len / pred_len, ov_len / gene_len)
            if recip > best_ov:
                best_ov, best_idx = recip, i
        if best_ov >= 0.8 and best_idx is not None:
            tp += 1
            matched[best_idx] = True
            if abs(pred["start"] - known[best_idx]["start"]) <= 3:
                correct_starts += 1
        else:
            fp += 1

    fn_ = sum(1 for m in matched if not m)
    n_pred, n_known = len(predicted), len(known)
    sens = tp / n_known if n_known else 0
    prec = tp / n_pred if n_pred else 0
    f1 = 2 * sens * prec / (sens + prec) if (sens + prec) > 0 else 0
    start_acc = correct_starts / tp if tp > 0 else 0

    return {
        "predicted": n_pred, "known": n_known,
        "tp": tp, "fp": fp, "fn": fn_,
        "sensitivity": sens, "precision": prec, "f1": f1,
        "start_accuracy": start_acc,
    }


def parse_gff_cds(path):
    genes = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            genes.append({
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
            })
    return sorted(genes, key=lambda g: g["start"])


def fmt(m: dict) -> str:
    return (f"Pred:{m['predicted']} Known:{m['known']} "
            f"TP:{m['tp']} FP:{m['fp']} FN:{m['fn']} "
            f"Sens:{m['sensitivity']*100:.1f}% Prec:{m['precision']*100:.1f}% "
            f"F1:{m['f1']:.3f} Start:{m['start_accuracy']*100:.1f}%")


def build_prompt(fasta_path: str, state: dict, cycle: int,
                 compile_error: str | None = None) -> str:
    p = f"""Write a COMPLETE standalone Rust program that annotates genes in a bacterial genome.

The program must:
1. Read FASTA file from command line arg: {fasta_path}
2. Find all genes using self-training hexamer models
3. Print JSON array to stdout: [{{"start":N,"end":M,"strand":"+"}}]
4. 1-based inclusive coordinates. Only std library.

Current best: src/bin/gene_finder.rs (read it and IMPROVE it)
Best F1: {state.get('best_f1', 0):.3f}
Target: F1 > 0.950 (Prodigal baseline)

"""
    if compile_error:
        p += f"PREVIOUS CODE FAILED TO COMPILE:\n{compile_error}\n\nFix the errors.\n\n"

    if cycle >= 1:
        ideas_path = WORK_DIR / "ideas.md"
        if ideas_path.exists():
            p += f"Ideas to incorporate: read {ideas_path}\n\n"

    reflect_path = WORK_DIR / f"reflection_{cycle-1}.txt"
    if reflect_path.exists():
        p += f"Previous reflection: read {reflect_path}\n\n"

    p += "Output the COMPLETE Rust source as JSON: {\"code\": \"use std...\"}\n"
    return p


def run_pipeline(
    fasta_path: str, gff_path: str,
    cross_fasta: str | None = None, cross_gff: str | None = None,
    n_cycles: int = 20,
    state_path: str = "annotator_work/rust_research_state.json",
):
    WORK_DIR.mkdir(exist_ok=True)

    # Ensure candidate bin exists in Cargo.toml
    cargo = Path("Cargo.toml").read_text()
    if "gene_finder_candidate" not in cargo:
        cargo += '\n[[bin]]\nname = "gene_finder_candidate"\npath = "src/bin/gene_finder_candidate.rs"\n'
        Path("Cargo.toml").write_text(cargo)

    try:
        state = json.loads(Path(state_path).read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        state = {"best_f1": 0.0, "best_cycle": -1, "metric_history": [], "current_cycle": 0}

    start_cycle = state.get("current_cycle", 0)
    fasta_abs = str(Path(fasta_path).resolve())
    gff_abs = str(Path(gff_path).resolve())

    print(f"═══ Rust Gene Finder Pipeline (cycles {start_cycle}→{start_cycle + n_cycles}) ═══")
    print(f"  Primary: {fasta_path}")
    if cross_fasta:
        print(f"  Cross-val: {cross_fasta}")
    print(f"  Best F1: {state.get('best_f1', 0):.3f}")
    print()

    for i in range(n_cycles):
        cycle = start_cycle + i
        print(f"━━━ Cycle {cycle} ━━━")
        t0 = time.time()

        # 1. GENERATE Rust code
        print("  [1/5] Generating Rust code...")
        prompt = build_prompt(fasta_abs, state, cycle)
        response = call_llm(prompt)
        code = extract_rust_code(response)
        if not code:
            print("  ✗ No Rust code from LLM")
            state["current_cycle"] = cycle + 1
            Path(state_path).write_text(json.dumps(state, indent=2))
            print()
            continue

        RUST_SRC.write_text(code)
        (WORK_DIR / f"rust_code_{cycle}.rs").write_text(code)
        print(f"  {len(code.splitlines())} lines of Rust")

        # 2. COMPILE
        print("  [2/5] Compiling...")
        ok, errors = compile_rust(RUST_SRC)
        if not ok:
            print(f"  ✗ Compile failed: {errors[:200]}")
            # Retry with error feedback
            prompt2 = build_prompt(fasta_abs, state, cycle, compile_error=errors)
            response2 = call_llm(prompt2)
            code2 = extract_rust_code(response2)
            if code2:
                RUST_SRC.write_text(code2)
                ok2, errors2 = compile_rust(RUST_SRC)
                if ok2:
                    code = code2
                    print("  ✓ Fixed on retry")
                else:
                    print(f"  ✗ Still fails: {errors2[:150]}")
                    state["current_cycle"] = cycle + 1
                    Path(state_path).write_text(json.dumps(state, indent=2))
                    print()
                    continue
            else:
                state["current_cycle"] = cycle + 1
                Path(state_path).write_text(json.dumps(state, indent=2))
                print()
                continue

        # 3. RUN + EVALUATE
        print("  [3/5] Running on primary genome...")
        try:
            stdout, elapsed = run_rust_binary(fasta_abs)
            metrics = evaluate_predictions(stdout, gff_abs)
            if not metrics:
                print("  ✗ Can't parse output")
                continue
            print(f"  {fmt(metrics)} ({elapsed:.1f}s)")
        except Exception as e:
            print(f"  ✗ Runtime error: {e}")
            state["current_cycle"] = cycle + 1
            Path(state_path).write_text(json.dumps(state, indent=2))
            print()
            continue

        # 4. CROSS-VALIDATE
        cross_metrics = None
        if cross_fasta and cross_gff:
            print("  [4/5] Cross-validating...")
            try:
                # Replace fasta path in source and recompile
                cross_abs = str(Path(cross_fasta).resolve())
                cross_gff_abs = str(Path(cross_gff).resolve())
                stdout2, elapsed2 = run_rust_binary(cross_abs)
                cross_metrics = evaluate_predictions(stdout2, cross_gff_abs)
                if cross_metrics:
                    print(f"  Cross: {fmt(cross_metrics)} ({elapsed2:.1f}s)")
                    if cross_metrics["f1"] < metrics["f1"] - 0.1:
                        print("  ⚠ Overfitting!")
            except Exception as e:
                print(f"  ✗ Cross-val error: {e}")
        else:
            print("  [4/5] Cross-val skipped")

        # 5. UPDATE
        is_best = metrics["f1"] > state.get("best_f1", 0)
        if is_best:
            state["best_f1"] = metrics["f1"]
            state["best_cycle"] = cycle
            BEST_RUST.write_text(code)
            Path("best_gene_finder.rs").write_text(code)
            print(f"  [5/5] ★ New best! F1={metrics['f1']:.3f}")
        else:
            print(f"  [5/5] Not best ({metrics['f1']:.3f} vs {state['best_f1']:.3f})")

        # Reflect (every 3rd cycle)
        if cycle % 3 == 0:
            print("  Reflecting...")
            hist = state.get("metric_history", [])[-10:]
            hist_str = "\n".join(f"  c{h['cycle']}: F1={h['f1']:.3f}" for h in hist)
            ref_prompt = f"Analyze gene finder results:\nCurrent: {fmt(metrics)}\nHistory:\n{hist_str}\nWhat to improve? Be concise."
            reflection = call_llm(ref_prompt, timeout=60) or ""
            (WORK_DIR / f"reflection_{cycle}.txt").write_text(reflection[:1000])
            print(f"  {reflection[:150]}")

        # Save state
        entry = {"cycle": cycle, "f1": metrics["f1"],
                 "sensitivity": metrics["sensitivity"], "precision": metrics["precision"],
                 "predicted": metrics["predicted"], "time_s": elapsed}
        if cross_metrics:
            entry["cross_f1"] = cross_metrics["f1"]
        state.setdefault("metric_history", []).append(entry)
        state["current_cycle"] = cycle + 1
        Path(state_path).write_text(json.dumps(state, indent=2))

        total_time = time.time() - t0
        print(f"  Cycle time: {total_time:.0f}s")
        print()

    print(f"═══ Done. Best F1: {state['best_f1']:.3f} (cycle {state.get('best_cycle', '?')}) ═══")


if __name__ == "__main__":
    fasta = sys.argv[1] if len(sys.argv) > 1 else "../data/ecoli_k12.fasta"
    gff = sys.argv[2] if len(sys.argv) > 2 else "../data/ecoli_k12.gff"
    cross_fasta = sys.argv[3] if len(sys.argv) > 3 else "../data/salmonella_lt2.fasta"
    cross_gff = sys.argv[4] if len(sys.argv) > 4 else "../data/salmonella_lt2.gff"
    n = int(sys.argv[5]) if len(sys.argv) > 5 else 20

    run_pipeline(fasta, gff, cross_fasta, cross_gff, n)
