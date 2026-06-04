#!/usr/bin/env python3
"""Compare gold_runner predictions (E_pred=QP, E_KS) to ARPES references.
Usage: validate.py <pred.csv> <arpes_reference.csv>"""
import csv, math, sys

pred = {}
with open(sys.argv[1]) as f:
    for r in csv.DictReader(f):
        pred.setdefault(int(r["point_id"]), []).append(
            (float(r["E_pred_eV"]), float(r["E_KS_eV"])))
exp = {}
with open(sys.argv[2]) as f:
    for r in csv.DictReader(f):
        exp[int(r["point_id"])] = float(r["E_expt_eV"])

ids = sorted(set(pred) & set(exp))
def rmse(a): return math.sqrt(sum(x * x for x in a) / len(a))
# nearest predicted band (lowest only here, but keep general)
qp = [min(p[0] for p in pred[i]) if False else
      min((p[0] for p in pred[i]), key=lambda v: abs(v - exp[i])) for i in ids]
ks = [min((p[1] for p in pred[i]), key=lambda v: abs(v - exp[i])) for i in ids]
dqp = [qp[j] - exp[ids[j]] for j in range(len(ids))]
dks = [ks[j] - exp[ids[j]] for j in range(len(ids))]
print(f"  matched {len(ids)} points")
print(f"  KS vs ARPES: RMSE={rmse(dks):.4f} eV  mean={sum(dks)/len(dks):+.4f}  max|.|={max(abs(x) for x in dks):.3f}")
print(f"  QP vs ARPES: RMSE={rmse(dqp):.4f} eV  mean={sum(dqp)/len(dqp):+.4f}  max|.|={max(abs(x) for x in dqp):.3f}")
print(f"  improvement: {rmse(dks)/rmse(dqp):.1f}x")
