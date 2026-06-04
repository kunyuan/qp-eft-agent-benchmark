#!/usr/bin/env python3
"""Hidden-set validator for the EFT-QP agent benchmark."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


PASS_RMSE_EV = 0.20
PARTIAL_RMSE_EV = 0.40
DEFAULT_MAX_BANDS_PER_POINT = {
    "Na": 1,
    "K": 1,
    "Mg": 2,
    "Al": 3,
}


def read_reference(path: Path) -> dict[int, float]:
    with path.open(newline="") as f:
        rows = csv.DictReader(f)
        return {int(row["point_id"]): float(row["E_expt_eV"]) for row in rows}


def read_element_config(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def max_bands_per_point(config: dict[str, object], element: str) -> int | None:
    validation = config.get("validation")
    if isinstance(validation, dict) and "max_bands_per_point" in validation:
        return int(validation["max_bands_per_point"])
    return DEFAULT_MAX_BANDS_PER_POINT.get(element)


def read_predictions(path: Path) -> dict[int, list[float]]:
    predictions: dict[int, list[float]] = {}
    with path.open(newline="") as f:
        rows = csv.DictReader(f)
        required = {"point_id", "E_pred_eV"}
        if rows.fieldnames is None or not required.issubset(set(rows.fieldnames)):
            raise ValueError(f"{path} must contain columns {sorted(required)}")
        for row in rows:
            value = (row.get("E_pred_eV") or "").strip()
            if not value:
                continue
            predictions.setdefault(int(row["point_id"]), []).append(float(value))
    return predictions


def validate_prediction_shape(
    element: str,
    predictions: dict[int, list[float]],
    max_bands: int | None,
) -> None:
    if max_bands is None:
        return
    for point_id, values in sorted(predictions.items()):
        if len(values) > max_bands:
            raise ValueError(
                f"too many bands for {element} point_id={point_id}: "
                f"got {len(values)}, allowed at most {max_bands}"
            )


def verdict(rmse: float | None) -> str:
    if rmse is None:
        return "NO_PREDICTION"
    if rmse < PASS_RMSE_EV:
        return "PASS"
    if rmse < PARTIAL_RMSE_EV:
        return "PARTIAL"
    return "FAIL"


def score_element(reference: dict[int, float], predictions: dict[int, list[float]]) -> dict[str, object]:
    scored = sorted(set(reference) & set(predictions))
    missing = sorted(set(reference) - set(predictions))
    if not scored:
        return {
            "n_scored": 0,
            "n_missing": len(missing),
            "rmse_eV": None,
            "mae_eV": None,
            "mean_signed_eV": None,
            "max_abs_error_eV": None,
            "verdict": "NO_PREDICTION",
        }

    residuals = []
    for point_id in scored:
        ref = reference[point_id]
        nearest = min(predictions[point_id], key=lambda value: abs(value - ref))
        residuals.append(nearest - ref)

    rmse = math.sqrt(sum(value * value for value in residuals) / len(residuals))
    mae = sum(abs(value) for value in residuals) / len(residuals)
    mean_signed = sum(residuals) / len(residuals)
    max_abs = max(abs(value) for value in residuals)
    return {
        "n_scored": len(scored),
        "n_missing": len(missing),
        "rmse_eV": round(rmse, 6),
        "mae_eV": round(mae, 6),
        "mean_signed_eV": round(mean_signed, 6),
        "max_abs_error_eV": round(max_abs, 6),
        "verdict": verdict(rmse),
    }


def run_submission(submission_dir: Path, hidden_dir: Path) -> dict[str, object]:
    runner = submission_dir / "run_qp.py"
    if not runner.exists():
        raise FileNotFoundError(f"submission must provide {runner}")

    per_element: dict[str, object] = {}
    with tempfile.TemporaryDirectory(prefix="qp_eft_validation_") as tmp:
        out_dir = Path(tmp)
        for element_dir in sorted(path for path in hidden_dir.iterdir() if path.is_dir()):
            element = element_dir.name
            config = read_element_config(element_dir / "element_config.json")
            out_csv = out_dir / f"{element}_qp_bands.csv"
            cmd = [
                sys.executable,
                str(runner),
                "--element-config",
                str(element_dir / "element_config.json"),
                "--grid",
                str(element_dir / "grid.csv"),
                "--out",
                str(out_csv),
            ]
            proc = subprocess.run(cmd, cwd=submission_dir, text=True, capture_output=True)
            if proc.returncode != 0:
                raise RuntimeError(
                    f"runner failed for {element} with exit code {proc.returncode}\n"
                    f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
                )
            reference = read_reference(element_dir / "arpes_reference.csv")
            predictions = read_predictions(out_csv)
            validate_prediction_shape(element, predictions, max_bands_per_point(config, element))
            per_element[element] = score_element(reference, predictions)

    rmses = [
        result["rmse_eV"]
        for result in per_element.values()
        if isinstance(result, dict) and result["rmse_eV"] is not None
    ]
    mean_rmse = sum(rmses) / len(rmses) if rmses else None
    return {
        "per_element": per_element,
        "overall": {
            "mean_rmse_eV": round(mean_rmse, 6) if mean_rmse is not None else None,
            "verdict": verdict(mean_rmse),
            "thresholds_eV": {
                "pass_below": PASS_RMSE_EV,
                "partial_below": PARTIAL_RMSE_EV,
            },
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--submission-dir", required=True, type=Path)
    parser.add_argument("--hidden-dir", type=Path, default=Path(__file__).resolve().parent / "hidden")
    parser.add_argument("--json", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        result = run_submission(args.submission_dir.resolve(), args.hidden_dir.resolve())
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
