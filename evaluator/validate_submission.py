#!/usr/bin/env python3
"""Hidden-set validator for the EFT-QP agent benchmark.

Scores a submitted `run_qp.py` against held-out ARPES, with three guards that
make "low RMSE" actually mean "did the frozen-core physics":

1. Flooding guard — predictions are capped to the gold occupied-band count per
   k-point. An agent that floods the energy window to make nearest-band matching
   trivial is REJECTED (this was the fatal exploit in the original scorer).
2. Nearest-band matching is then safe: with n_pred ~= n_occ (1-3 bands), there is
   nothing to flood.
3. KS-baseline gate — calibrated thresholds (PASS<0.30 eV) sit below what bare
   uncorrected KS scores (~0.41-0.61 eV), so a no-physics submission FAILs. The
   report also states what KS-only would have scored, for audit.

Thresholds calibrated from the gold reference vs real ARPES (see DESIGN.md §5b).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

PASS_RMSE_EV = 0.30
PARTIAL_RMSE_EV = 0.40
# fraction of points allowed to exceed the gold band count before we call it flooding
FLOOD_POINT_FRAC = 0.20
# total predicted bands may not exceed this multiple of the gold occupied total
FLOOD_TOTAL_RATIO = 1.5


def read_reference(path: Path) -> dict[int, float]:
    with path.open(newline="") as f:
        return {int(r["point_id"]): float(r["E_expt_eV"]) for r in csv.DictReader(f)}


def copy_runner_inputs(element_dir: Path, dst: Path) -> None:
    """Copy the hidden runner inputs EXCEPT the held-out ARPES (the answer key),
    so a malicious run_qp.py cannot read the reference energies sitting next to
    the config it is handed. Keeps `element` and everything else the runner needs."""
    dst.mkdir(parents=True, exist_ok=True)
    for src in sorted(element_dir.iterdir()):
        if src.name == "arpes_reference.csv":
            continue
        target = dst / src.name
        if src.is_dir():
            shutil.copytree(src, target)
        else:
            shutil.copy2(src, target)


def read_predictions(path: Path) -> dict[int, list[float]]:
    out: dict[int, list[float]] = {}
    with path.open(newline="") as f:
        rows = csv.DictReader(f)
        required = {"point_id", "E_pred_eV"}
        if rows.fieldnames is None or not required.issubset(set(rows.fieldnames)):
            raise ValueError(f"{path} must contain columns {sorted(required)}")
        for r in rows:
            v = (r.get("E_pred_eV") or "").strip()
            if v:
                out.setdefault(int(r["point_id"]), []).append(float(v))
    return out


def read_gold(path: Path) -> tuple[dict[int, list[float]], dict[int, list[float]]]:
    """gold/<El>_gold.csv -> (qp_bands_per_point, ks_bands_per_point)."""
    qp: dict[int, list[float]] = {}
    ks: dict[int, list[float]] = {}
    with path.open(newline="") as f:
        for r in csv.DictReader(f):
            pid = int(r["point_id"])
            qp.setdefault(pid, []).append(float(r["E_pred_eV"]))
            ks.setdefault(pid, []).append(float(r["E_KS_eV"]))
    return qp, ks


def verdict(rmse: float | None) -> str:
    if rmse is None:
        return "NO_PREDICTION"
    if rmse < PASS_RMSE_EV:
        return "PASS"
    if rmse < PARTIAL_RMSE_EV:
        return "PARTIAL"
    return "FAIL"


def _rmse(residuals: list[float]) -> float:
    return math.sqrt(sum(x * x for x in residuals) / len(residuals))


def _nearest_rmse(reference: dict[int, float], bands: dict[int, list[float]]) -> float | None:
    res = []
    for pid, ref in reference.items():
        if pid in bands and bands[pid]:
            res.append(min(bands[pid], key=lambda v: abs(v - ref)) - ref)
    return _rmse(res) if res else None


def score_element(
    reference: dict[int, float],
    predictions: dict[int, list[float]],
    gold_qp: dict[int, list[float]],
    gold_ks: dict[int, list[float]],
) -> dict[str, object]:
    n_occ = {pid: len(b) for pid, b in gold_qp.items()}

    # --- flooding guard ---
    over = [pid for pid, b in predictions.items() if len(b) > n_occ.get(pid, 1) + 1]
    total_pred = sum(len(b) for b in predictions.values())
    total_gold = sum(n_occ.values()) or 1
    flooding = (
        len(over) > FLOOD_POINT_FRAC * max(len(predictions), 1)
        or total_pred > FLOOD_TOTAL_RATIO * total_gold
    )

    scored = sorted(set(reference) & set(predictions))
    missing = sorted(set(reference) - set(predictions))
    if flooding:
        return {
            "verdict": "REJECTED_FLOODING",
            "rmse_eV": None,
            "n_scored": len(scored),
            "n_missing": len(missing),
            "n_points_over_band_cap": len(over),
            "total_predicted_bands": total_pred,
            "total_gold_bands": total_gold,
            "note": "predictions exceed gold occupied-band count; cannot pass",
        }
    # shape guard: every reference point must be predicted with exactly the gold
    # band count (occupied + first unoccupied). Closes "sparse" (drop hard points)
    # and "under-band" cheats; honest agents on the pinned setup match exactly.
    band_count_mismatches = {
        pid: {"expected": n_occ.get(pid, 0), "got": len(predictions.get(pid, []))}
        for pid in reference
        if pid in predictions and len(predictions[pid]) != n_occ.get(pid, 0)
    }
    if missing or band_count_mismatches:
        return {
            "verdict": "INVALID_SHAPE",
            "rmse_eV": None,
            "n_scored": len(scored),
            "n_missing": len(missing),
            "n_band_count_mismatches": len(band_count_mismatches),
            "note": "every reference point must be predicted with exactly the gold band count",
        }
    if not scored:
        return {"verdict": "NO_PREDICTION", "rmse_eV": None,
                "n_scored": 0, "n_missing": len(missing)}

    residuals = [min(predictions[pid], key=lambda v: abs(v - reference[pid])) - reference[pid]
                 for pid in scored]
    rmse = _rmse(residuals)
    return {
        "verdict": verdict(rmse),
        "rmse_eV": round(rmse, 6),
        "mae_eV": round(sum(abs(x) for x in residuals) / len(residuals), 6),
        "mean_signed_eV": round(sum(residuals) / len(residuals), 6),
        "max_abs_error_eV": round(max(abs(x) for x in residuals), 6),
        "n_scored": len(scored),
        "n_missing": len(missing),
        "n_points_over_band_cap": len(over),
        # audit: what bare KS would have scored (the gate the agent must beat)
        "ks_baseline_rmse_eV": round(r, 6) if (r := _nearest_rmse(reference, gold_ks)) else None,
    }


def run_submission(submission_dir: Path, hidden_dir: Path, gold_dir: Path) -> dict[str, object]:
    runner = submission_dir / "run_qp.py"
    if not runner.exists():
        raise FileNotFoundError(f"submission must provide {runner}")

    per_element: dict[str, object] = {}
    with tempfile.TemporaryDirectory(prefix="qp_eft_validation_") as tmp:
        out_dir = Path(tmp)
        for element_dir in sorted(p for p in hidden_dir.iterdir() if p.is_dir()):
            el = element_dir.name
            # hand the runner a sanitized input dir with NO arpes_reference.csv
            runner_input = out_dir / "inputs" / el
            copy_runner_inputs(element_dir, runner_input)
            out_csv = out_dir / f"{el}_qp_bands.csv"
            cmd = [sys.executable, str(runner),
                   "--element-config", str(runner_input / "element_config.json"),
                   "--grid", str(runner_input / "grid.csv"),
                   "--out", str(out_csv)]
            proc = subprocess.run(cmd, cwd=submission_dir, text=True, capture_output=True)
            if proc.returncode != 0:
                raise RuntimeError(
                    f"runner failed for {el} (exit {proc.returncode})\n"
                    f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}")
            reference = read_reference(element_dir / "arpes_reference.csv")
            predictions = read_predictions(out_csv)
            gold_qp, gold_ks = read_gold(gold_dir / f"{el}_gold.csv")
            per_element[el] = score_element(reference, predictions, gold_qp, gold_ks)

    rmses = [r["rmse_eV"] for r in per_element.values()
             if isinstance(r, dict) and r.get("rmse_eV") is not None]
    mean_rmse = sum(rmses) / len(rmses) if rmses else None

    def _has(v):
        return any(isinstance(r, dict) and r.get("verdict") == v for r in per_element.values())
    # any per-element disqualifier sinks the whole submission
    if _has("REJECTED_FLOODING"):
        overall = "REJECTED_FLOODING"
    elif _has("INVALID_SHAPE"):
        overall = "INVALID_SHAPE"
    elif _has("NO_PREDICTION"):
        overall = "NO_PREDICTION"
    else:
        overall = verdict(mean_rmse)
    return {
        "per_element": per_element,
        "overall": {
            "mean_rmse_eV": round(mean_rmse, 6) if mean_rmse is not None else None,
            "verdict": overall,
            "thresholds_eV": {"pass_below": PASS_RMSE_EV, "partial_below": PARTIAL_RMSE_EV},
        },
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--submission-dir", required=True, type=Path)
    p.add_argument("--level", type=int, choices=(1, 2, 3), default=2,
                   help="difficulty level; selects hidden/L<level> (default 2)")
    p.add_argument("--hidden-dir", type=Path, default=None,
                   help="override; defaults to evaluator/hidden/L<level>")
    p.add_argument("--gold-dir", type=Path, default=Path(__file__).resolve().parent / "gold")
    p.add_argument("--json", type=Path)
    args = p.parse_args()
    if args.hidden_dir is None:
        args.hidden_dir = Path(__file__).resolve().parent / "hidden" / f"L{args.level}"
    return args


def main() -> int:
    args = parse_args()
    try:
        result = run_submission(args.submission_dir.resolve(),
                                args.hidden_dir.resolve(), args.gold_dir.resolve())
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        return 2
    text = json.dumps(result, indent=2, sort_keys=True)
    print(text)
    if args.json:
        args.json.write_text(text + "\n")
    return 0 if result["overall"]["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
