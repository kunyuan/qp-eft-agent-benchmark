#!/usr/bin/env python3
"""Per-band before/after diagnostic: bare DFT (KS) vs frozen-core QP, each
compared to experiment under the SAME band assignment.

Shows where the correction acts: deep bands (large KS overbinding) are pulled
into agreement, while bands near E_F (where z_core*(E_KS-E_F) -> 0) are left
essentially unchanged. A coherent per-band match also confirms the grid's
k-mapping is sound for multi-band elements (the experiment tracks real bands,
not a loose nearest-of-many artifact).

Usage:
  python validate_bands.py <El> [El ...]      # uses evaluator/gold + source_data
  python validate_bands.py --gold a.csv --arpes b.csv
"""
from __future__ import annotations
import argparse, csv, math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(gold_csv: Path, arpes_csv: Path):
    pairs = defaultdict(list)  # point_id -> [(QP, KS), ...] (row order = band order)
    for r in csv.DictReader(open(gold_csv)):
        pairs[int(r["point_id"])].append((float(r["E_pred_eV"]), float(r["E_KS_eV"])))
    arpes = {int(r["point_id"]): float(r["E_expt_eV"])
             for r in csv.DictReader(open(arpes_csv))}
    return pairs, arpes


def rmse(a):
    return math.sqrt(sum(x * x for x in a) / len(a)) if a else float("nan")


def diagnose(name: str, gold_csv: Path, arpes_csv: Path):
    pairs, arpes = load(gold_csv, arpes_csv)
    byband = defaultdict(lambda: {"ks": [], "qp": []})
    for pid, e in arpes.items():
        pr = sorted(pairs[pid])                       # deep -> shallow by QP
        qp, ks = min(pr, key=lambda p: abs(p[0] - e))  # match by QP; KS = same band
        bi = pr.index((qp, ks)) + 1
        byband[bi]["ks"].append(ks - e)
        byband[bi]["qp"].append(qp - e)

    print(f"\n=== {name}: deviation from experiment (eV), same band assignment ===")
    print(f"{'band':>5} {'n':>4} | {'DFT RMSE':>9} {'mean':>7} | {'QP RMSE':>8} {'mean':>7} | improve")
    allks, allqp = [], []
    for bi in sorted(byband):
        ks, qp = byband[bi]["ks"], byband[bi]["qp"]
        allks += ks
        allqp += qp
        imp = rmse(ks) / rmse(qp) if rmse(qp) else float("inf")
        print(f"{bi:5d} {len(ks):4d} | {rmse(ks):9.3f} {sum(ks)/len(ks):+7.3f} | "
              f"{rmse(qp):8.3f} {sum(qp)/len(qp):+7.3f} | {imp:.1f}x")
    imp = rmse(allks) / rmse(allqp) if rmse(allqp) else float("inf")
    print(f"{'ALL':>5} {len(allks):4d} | {rmse(allks):9.3f} {sum(allks)/len(allks):+7.3f} | "
          f"{rmse(allqp):8.3f} {sum(allqp)/len(allqp):+7.3f} | {imp:.1f}x")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("elements", nargs="*", help="element symbols (Na K Al Mg ...)")
    p.add_argument("--gold", type=Path)
    p.add_argument("--arpes", type=Path)
    a = p.parse_args()
    if a.gold and a.arpes:
        diagnose(a.gold.stem, a.gold, a.arpes)
        return
    for el in (a.elements or ["Na", "K", "Mg", "Al"]):
        diagnose(el, ROOT / "evaluator" / "gold" / f"{el}_gold.csv",
                 ROOT / "evaluator" / "source_data" / el / "arpes_reference.csv")


if __name__ == "__main__":
    main()
