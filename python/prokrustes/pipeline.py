"""Multi-stage research pipeline: assemble → test → cross-validate → reflect → review → update.

This is the main loop that replaces the simple generate→run→check cycle.
Each stage is a separate LLM call with a focused role.
"""

import json
import os
import sys
import subprocess
import time
from pathlib import Path

WORK_DIR = Path("./annotator_work")
LIB_DIR = Path("./gene_finder_lib")


def call_llm(prompt: str, timeout: int = 300) -> str | None:
    """Call claude with tools and permissions."""
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
    except (subprocess.TimeoutExpired, json.JSONDecodeError, Exception) as e:
        print(f"  LLM error: {e}")
        return None


def extract_code(response: str) -> str | None:
    """Extract Python code from LLM response."""
    if not response:
        return None
    # Try JSON
    for marker in ["```json", "```python"]:
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
    # Try raw JSON
    try:
        return json.loads(response).get("code")
    except (json.JSONDecodeError, AttributeError):
        pass
    # Find def/import
    for marker in ["import ", "def ", "#!/"]:
        if marker in response:
            return response[response.index(marker):]
    return None


def run_script(code: str, timeout: int = 120) -> tuple[str, str, int]:
    """Run Python script, return (stdout, stderr, returncode)."""
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
        f.write(code)
        tmp = f.name
    try:
        result = subprocess.run(
            ["python3", tmp], capture_output=True, text=True, timeout=timeout,
        )
        return result.stdout, result.stderr, result.returncode
    except subprocess.TimeoutExpired:
        return "", "TIMEOUT", 1
    finally:
        os.unlink(tmp)


def evaluate_predictions(pred_json: str, gff_path: str) -> dict | None:
    """Compare predictions with known genes. Returns metrics dict."""
    try:
        predicted = json.loads(pred_json)
    except json.JSONDecodeError:
        # Try finding JSON array in output
        if '[' in pred_json:
            try:
                predicted = json.loads(pred_json[pred_json.index('['):pred_json.rindex(']')+1])
            except (json.JSONDecodeError, ValueError):
                return None
        else:
            return None

    known = parse_gff_cds(gff_path)

    tp, fp, fn_ = 0, 0, 0
    correct_starts = 0
    matched = [False] * len(known)

    for pred in predicted:
        best_ov, best_idx = 0.0, None
        for i, gene in enumerate(known):
            if matched[i]:
                continue
            if pred.get("strand", "+") != gene["strand"]:
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
    n_pred = len(predicted)
    n_known = len(known)
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


def format_metrics(m: dict) -> str:
    return (f"Predicted: {m['predicted']} | Known: {m['known']} | "
            f"TP: {m['tp']} FP: {m['fp']} FN: {m['fn']} | "
            f"Sens: {m['sensitivity']*100:.1f}% Prec: {m['precision']*100:.1f}% "
            f"F1: {m['f1']:.3f} | Start acc: {m['start_accuracy']*100:.1f}%")


# ═══════════════════════════════════════════════════════════
# Pipeline stages
# ═══════════════════════════════════════════════════════════


def stage_assemble(cycle: int, state: dict, fasta_path: str) -> str | None:
    """Stage 1: LLM assembles gene finder from library parts or creates new."""
    lib_desc = ""
    for name in sorted(os.listdir(LIB_DIR)):
        if name.endswith(".py") and name != "__init__.py" and name != "pipeline.py" and name != "registry.py":
            lib_desc += f"\n  {LIB_DIR / name}"

    prompt = f"""You are building a bacterial gene annotator.

FASTA file: {fasta_path}

Available library parts (read them):{lib_desc}

Previous best code: {WORK_DIR / 'previous_code.py'}
Previous best F1: {state.get('best_f1', 0):.3f}
Target: F1 > 0.950 (Prodigal baseline)

IMPORTANT: Read the library parts first. Use them as building blocks.
If you need functionality not in the library, write it.
Output: complete standalone Python script as JSON {{"code": "..."}}
The script must read the FASTA file and print JSON gene array to stdout."""

    response = call_llm(prompt)
    if not response:
        return None
    return extract_code(response)


def stage_test(code: str, fasta_path: str, gff_path: str) -> dict | None:
    """Stage 2: Run on primary genome, compute metrics."""
    stdout, stderr, rc = run_script(code)
    if rc != 0:
        print(f"  ✗ Script error: {stderr[:300]}")
        return None
    metrics = evaluate_predictions(stdout, gff_path)
    if metrics is None:
        print(f"  ✗ Can't parse output: {stdout[:200]}")
    return metrics


def stage_cross_validate(code: str, fasta_path: str, gff_path: str,
                          original_fasta: str) -> dict | None:
    """Stage 3: Run on secondary genome WITHOUT retraining.

    Replace FASTA path in code and run. This tests generalization.
    """
    # Replace the hardcoded fasta path
    cross_code = code.replace(original_fasta, fasta_path)
    if cross_code == code:
        # Path not found as string — try common patterns
        for pattern in ['FASTA_PATH = "', "FASTA_PATH = '", 'fasta_path = "', "path = '"]:
            if pattern in code:
                old_line = code[code.index(pattern):code.index('\n', code.index(pattern))]
                new_line = f'{pattern.split("=")[0].strip()} = "{fasta_path}"'
                cross_code = code.replace(old_line, new_line)
                break

    stdout, stderr, rc = run_script(cross_code)
    if rc != 0:
        print(f"  ✗ Cross-validation error: {stderr[:200]}")
        return None
    return evaluate_predictions(stdout, gff_path)


def stage_reflect(cycle: int, state: dict, metrics: dict,
                   cross_metrics: dict | None) -> str:
    """Stage 4: LLM reflects on results — what worked, what didn't, why."""
    history = state.get("metric_history", [])
    history_str = "\n".join(
        f"  Cycle {h['cycle']}: F1={h['f1']:.3f} {h.get('action', '')}"
        for h in history[-10:]
    )

    cross_str = format_metrics(cross_metrics) if cross_metrics else "Not available"

    prompt = f"""You are a research scientist reflecting on gene finding results.

## Current results (cycle {cycle})
{format_metrics(metrics)}

## Cross-validation (Salmonella)
{cross_str}

## History (last 10 cycles)
{history_str}

## Best so far: F1={state.get('best_f1', 0):.3f}
## Target: F1 > 0.950 (Prodigal)

Analyze:
1. What improved compared to previous cycles?
2. What got worse?
3. Why? Look at TP/FP/FN changes.
4. If cross-validation F1 << primary F1, the model is overfitting.
5. What specific changes would help most?

Be concise and specific. No code — just analysis."""

    response = call_llm(prompt, timeout=60)
    return response or "No reflection available"


def stage_code_review(code: str) -> str:
    """Stage 5: LLM reviews code for bugs."""
    # Only send first/last parts to keep prompt small
    lines = code.split('\n')
    if len(lines) > 200:
        code_preview = '\n'.join(lines[:100]) + '\n...\n' + '\n'.join(lines[-100:])
    else:
        code_preview = code

    prompt = f"""Review this bacterial gene finder for bugs. Focus on:
1. Off-by-one errors in coordinates (1-based vs 0-based)
2. Strand handling errors (reverse complement)
3. Division by zero
4. Edge cases: genes at genome start/end, very short/long genes
5. Logic errors in scoring/filtering

Code:
```python
{code_preview}
```

List ONLY actual bugs (not style issues). For each bug: line description + fix.
If no bugs found, say "No bugs found."
Be concise."""

    response = call_llm(prompt, timeout=60)
    return response or "No review available"


def stage_update_library(code: str, metrics: dict, state: dict):
    """Stage 6: If new code has useful parts, extract to library."""
    # Only update if it's the new best
    if metrics["f1"] <= state.get("best_f1", 0):
        return

    # Save as best
    (WORK_DIR / "previous_code.py").write_text(code)
    Path("best_annotator.py").write_text(code)


# ═══════════════════════════════════════════════════════════
# Main pipeline
# ═══════════════════════════════════════════════════════════


def run_pipeline(
    fasta_path: str,
    gff_path: str,
    cross_fasta: str | None = None,
    cross_gff: str | None = None,
    n_cycles: int = 20,
    state_path: str = "annotator_work/research_state.json",
):
    """Run the full multi-stage pipeline."""
    WORK_DIR.mkdir(exist_ok=True)

    # Load state
    try:
        state = json.loads(Path(state_path).read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        state = {"best_f1": 0.0, "best_cycle": -1, "metric_history": [], "current_cycle": 0}

    start_cycle = state.get("current_cycle", 0)
    fasta_abs = str(Path(fasta_path).resolve())

    print(f"═══ Gene Finder Pipeline (cycles {start_cycle}→{start_cycle + n_cycles}) ═══")
    print(f"  Primary: {fasta_path}")
    if cross_fasta:
        print(f"  Cross-val: {cross_fasta}")
    print(f"  Best F1: {state.get('best_f1', 0):.3f}")
    print()

    for i in range(n_cycles):
        cycle = start_cycle + i
        print(f"━━━ Cycle {cycle} ━━━")

        # 1. ASSEMBLE
        print("  [1/6] Assembling...")
        code = stage_assemble(cycle, state, fasta_abs)
        if not code:
            print("  ✗ No code generated")
            state["current_cycle"] = cycle + 1
            Path(state_path).write_text(json.dumps(state, indent=2))
            continue

        (WORK_DIR / f"code_{cycle}.py").write_text(code)

        # 2. TEST
        print("  [2/6] Testing on primary genome...")
        metrics = stage_test(code, fasta_abs, gff_path)
        if not metrics:
            print("  ✗ Test failed")
            state["current_cycle"] = cycle + 1
            Path(state_path).write_text(json.dumps(state, indent=2))
            continue

        print(f"  {format_metrics(metrics)}")

        # 3. CROSS-VALIDATE
        cross_metrics = None
        if cross_fasta and cross_gff:
            print("  [3/6] Cross-validating on Salmonella...")
            cross_metrics = stage_cross_validate(
                code, str(Path(cross_fasta).resolve()),
                cross_gff, fasta_abs,
            )
            if cross_metrics:
                print(f"  Cross: {format_metrics(cross_metrics)}")
                if cross_metrics["f1"] < metrics["f1"] - 0.1:
                    print("  ⚠ Overfitting detected! Cross-val F1 much lower.")
        else:
            print("  [3/6] Cross-validation skipped (no secondary genome)")

        # 4. REFLECT
        print("  [4/6] Reflecting...")
        reflection = stage_reflect(cycle, state, metrics, cross_metrics)
        reflection_short = (reflection or "")[:200]
        print(f"  {reflection_short}")
        (WORK_DIR / f"reflection_{cycle}.txt").write_text(reflection or "")

        # 5. CODE REVIEW (every 3rd cycle to save LLM calls)
        if cycle % 3 == 0:
            print("  [5/6] Code review...")
            review = stage_code_review(code)
            review_short = (review or "")[:200]
            print(f"  {review_short}")
            (WORK_DIR / f"review_{cycle}.txt").write_text(review or "")
        else:
            print("  [5/6] Code review skipped (every 3rd cycle)")

        # 6. UPDATE
        is_best = metrics["f1"] > state.get("best_f1", 0)
        if is_best:
            stage_update_library(code, metrics, state)
            state["best_f1"] = metrics["f1"]
            state["best_cycle"] = cycle
            print(f"  [6/6] ★ New best! F1={metrics['f1']:.3f}")
        else:
            print(f"  [6/6] Not best (current: {metrics['f1']:.3f}, best: {state['best_f1']:.3f})")

        # Update state
        state["current_cycle"] = cycle + 1
        history_entry = {
            "cycle": cycle, "f1": metrics["f1"],
            "sensitivity": metrics["sensitivity"],
            "precision": metrics["precision"],
            "predicted": metrics["predicted"],
            "action": "best" if is_best else "try",
        }
        if cross_metrics:
            history_entry["cross_f1"] = cross_metrics["f1"]
        state.setdefault("metric_history", []).append(history_entry)
        Path(state_path).write_text(json.dumps(state, indent=2))
        print()

    print(f"═══ Done. Best F1: {state['best_f1']:.3f} (cycle {state.get('best_cycle', '?')}) ═══")


if __name__ == "__main__":
    fasta = sys.argv[1] if len(sys.argv) > 1 else "../data/ecoli_k12.fasta"
    gff = sys.argv[2] if len(sys.argv) > 2 else "../data/ecoli_k12.gff"
    cross_fasta = sys.argv[3] if len(sys.argv) > 3 else "../data/salmonella_lt2.fasta"
    cross_gff = sys.argv[4] if len(sys.argv) > 4 else "../data/salmonella_lt2.gff"
    n = int(sys.argv[5]) if len(sys.argv) > 5 else 20

    run_pipeline(fasta, gff, cross_fasta, cross_gff, n)
