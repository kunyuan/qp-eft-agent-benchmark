# L3 probe — agent "r5": analysis & reproduction

One of **6 parallel, no-web, fresh agents** asked to solve **Level 3** (derive the
frozen-core quasiparticle correction from the physical setup alone — no formula
given) and implement it as `run_qp.py`. This report analyses r5 (the strongest of
the six) and gives exact steps to reproduce its result.

Solution preserved here: [`solution/run_qp.py`](solution/run_qp.py),
[`solution/compute_bands.jl`](solution/compute_bands.jl),
[`solution/method.md`](solution/method.md) (r5's own derivation, verbatim).
Outputs it produced: [`outputs/`](outputs/).

---

## TL;DR

- r5 produced a **genuine, independent derivation** (a core-polarizability /
  closure-variance pole rescaling) — **not** the gold's formula, and different from
  all 5 other agents. Evidence it derived rather than recalled: 6 agents → 6
  distinct forms.
- On the **visible** dev elements (Na, Al) it looks **gold-quality** (Na 0.077 eV).
- On the **hidden** held-out metals it **fails**: K is **over-narrowed**, Mg is
  **under-narrowed**, and Mg also trips a band-count bug → overall **INVALID_SHAPE
  (not pass)**.
- It is the **best of the six on K** (lowest K RMSE) and structurally the closest to
  the gold, but it is **not** the most accurate overall (its Mg has no valid number).
- **Root cause (physics):** the *energy structure* (`1/ΔE_c²` pole rescaling) is
  right; the *coupling / form factor* is wrong — a real-space variance overlapped
  with the pseudo valence density, instead of the gold's momentum-space form factor.
  It over-couples diffuse cores (K) and under-couples tight cores (Mg).
- **Root cause (engineering):** the real-space nearest-atom binning + `eigs < eF`
  band-emit logic mis-counts bands at one k-point on the 2-atom hcp Mg cell.

---

## 1. What r5 derived (its physics)

Full text in [`solution/method.md`](solution/method.md). The chain:

1. The large-core PSP keeps the **static** core Hartree field `V_H_c`; what is
   missing is the **dynamical** core response (a valence electron virtually excites
   a core s-shell `c→c'` at cost `ΔE_c`) — a 2nd-order self-energy.
2. **Closure (Unsöld)** over the core manifold at the single scale `ΔE_c` turns the
   excited-state sum into the **variance of the core Coulomb field**:
   `F_c(r) = ⟨c|1/|r−r'|²|c⟩ − V_H_c(r)²`.
3. By the Galilean-invariance hint, the net effect is a **pole rescaling**:
   `E_QP−E_F = Z_nk·(ε_nk−E_F)`, `Z_nk = 1/(1+λ_nk)`.
4. **Final closed form:**

   ```
   λ_nk = (1/π) · Σ_c N_c · (Ω/n_atom)·⟨ψ_nk|F_c|ψ_nk⟩ / ΔE_c²,   N_c = 2
   ```

   with `⟨ψ_nk|F_c|ψ_nk⟩` evaluated by transforming `c_nk(G)` to real space and
   **binning `|ψ_nk(r)|²` by distance to the nearest atom**.

r5's own caveat: the `1/π` and `Ω/n_atom` prefactors are "fixed at the **joint
optimum** for both public elements" — i.e. calibrated to Na/Al.

---

## 2. Results

### Occupied bandwidth (Γ-point depth below E_F, eV) — r5 vs references

| element | KS (bare) | **r5 QP** | gold QP | eDMFT | ARPES | r5 Z | correct Z |
|---------|-----------|-----------|---------|-------|-------|------|-----------|
| Na (dev)  | 3.266  | **2.676** | 2.616  | 2.84 | 2.65–2.78 | 0.819 | 0.801 |
| Al (dev)  | 11.179 | **10.683**| 10.704 | —    | ~10.6     | 0.956 | 0.957 |
| **K** (hidden)  | 2.270  | **1.226** | 1.492 | 1.42 | 1.60 | **0.540** | 0.657 |
| **Mg** (hidden) | 6.921  | **6.616** | 6.310 | 6.18 | 6.15 | **0.956** | 0.912 |

- Na/Al: excellent (Na sits inside the ARPES window).
- **K: over-narrowed** (Z 0.540 ≪ 0.657) — band 0.37 eV too shallow vs ARPES.
- **Mg: under-narrowed** (Z 0.956 ≫ 0.912) — band 0.47 eV too deep vs ARPES.

### Scored through the evaluator (`--level 3`, tight bars K<0.17, Mg<0.21)

| element | RMSE vs ARPES | verdict |
|---------|---------------|---------|
| K  | 0.260 | PARTIAL (over the 0.17 bar) |
| Mg | — | **INVALID_SHAPE** (1 missing point, 1 band-count mismatch) |
| overall | — | **INVALID_SHAPE** → not pass |

---

## 3. Where it goes wrong

### (A) Physics: the coupling/form factor, not the energy denominator

Write the correction strength `λ = 1/Z − 1` and compare to the gold and to a
hypothetical *element-universal* coupling (`λ ∝ 1/ΔE_c²`, constant fit to Na):

| el | ΔE_c (Ha) | gold λ | **r5 λ** | universal-W λ |
|----|-----------|--------|----------|---------------|
| Na | 2.70 | 0.248 | 0.220 | 0.248 |
| Al | 5.69 | 0.045 | 0.046 | 0.056 |
| K  | 1.72 | 0.522 | **0.852** ← over | 0.612 |
| Mg | 4.07 | 0.096 | **0.046** ← under | 0.110 |

Reading: the `1/ΔE_c²` structure (universal-W column) already lands K/Mg **near**
the gold. r5's **actual** λ is *further* from the gold than that simple trend — so
the error lives in r5's **per-element coupling** `W_c = ⟨ψ|F_c|ψ⟩`, not in `ΔE_c`.

Concretely, r5 uses the **real-space** Coulomb-field variance
`F_c(r) = ⟨1/|r−r'|²⟩_c − V_H_c²` overlapped with the **pseudo** valence density
`|ψ_nk|²`. The gold uses the **momentum-space form factor**
`f_c(K) = √(4π)/K ∫ u_c(r)[V_H_c(r) − J_c] sin(Kr) dr`, evaluated at `K=|k+G|` and
contracted with `c_nk(G)`. These are different operators; the real-space-variance
approximation **over-estimates** the coupling for K's diffuse 3s core and
**under-estimates** it for Mg's tight 2s core. Calibrating the `1/π` prefactor to
Na/Al then **locks in** the visible-element fit and pushes the residual onto K/Mg.

> Irony: r5's own insight "`W_c` is near-universal (~0.8 Ha)" was *better* than its
> detailed computation. Had it trusted `λ ∝ Σ 1/ΔE_c²` (universal-W column), K and
> Mg would have come out much closer. Its mistake was computing the real-space
> overlap per element — exactly the miscalibrated ingredient.

### (B) Engineering: Mg band-count (INVALID_SHAPE)

A separate bug. Band emission in `run_qp.py` is `occ = (eigs < eF)`, emit
`occupied + first unoccupied`. On the **hcp 2-atom Mg** cell, the real-space
nearest-atom binning + this `eigs<eF` count produce the wrong band count at one
k-point (`n_missing=1, n_band_count_mismatches=1`) → the whole element is rejected.
Even with perfect physics Mg would have been thrown out on shape. (K already fails,
so fixing Mg would not change the overall verdict.)

---

## 4. Standing among the 6 attempts (hidden K/Mg)

| run | K | Mg | overall | derivation |
|-----|------|------|---------|-----------|
| r1 | 0.325 | 0.244 | FAIL | `1/(1+Σ(J_c/ΔE_c)²)` |
| r2 | 0.420 | **0.212** | FAIL | state-dep `⟨r²⟩⟨V_Hc⟩_nk/ΔE²` |
| r3 | 0.768 | 0.228 | FAIL | `1−Σ(2/3)⟨r²⟩/ΔE` |
| r4 | crash | crash | CRASHED | state-dep `f_l0(nk)` |
| **r5** | **0.260** | INVALID | INVALID_SHAPE | closure variance `F_c`, `1/π` |
| r6 | crash | crash | CRASHED | `1−Σ J_c/ΔE_c²` |

- **Best K:** r5 (0.260). **Best Mg:** r2 (0.212, within 0.002 of its bar).
- **Most complete:** r1 (both valid). **Solve rate: 0/6.**
- No single agent got *both* K and Mg right — the L3 difficulty is fitting one
  parameter-free expression to two opposite extremes (diffuse vs tight core).

---

## 5. Reproduction

From the repo root, with Julia + the pinned project (`environment/`) available:

```bash
export JULIA_PROJECT="$PWD/environment"      # r5's pinned DFTK 0.7.25 env
export OPENBLAS_NUM_THREADS=2                 # keep it polite on shared machines

# (a) run r5's solver on any single element ------------------------------------
#   public dev cases:  agent_packet/levels/L3/{Na,Al}/
#   hidden cases:      evaluator/hidden/L3/{K,Mg}/   (maintainer-only)
python3 experiments/L3_r5/solution/run_qp.py \
  --element-config evaluator/hidden/L3/K/element_config.json \
  --grid          evaluator/hidden/L3/K/grid.csv \
  --out           /tmp/r5_K.csv
# QP Γ-depth = -min(E_pred_eV) over the t=0 point. Expected ~1.23 eV (over-narrow).

# (b) score on the full hidden set through the evaluator ------------------------
python3 evaluator/validate_submission.py \
  --submission-dir experiments/L3_r5/solution --level 3 --json /tmp/r5_score.json
# Expected: K PARTIAL ~0.260, Mg INVALID_SHAPE, overall INVALID_SHAPE.
```

Expected QP Γ-depths (reproduces [`outputs/`](outputs/)): Na 2.676, Al 10.683,
K 1.226, Mg 6.616 eV.

Notes:
- `run_qp.py` reads `core_model.json` + `atomic_core_*.csv` **relative to the
  config path**, so any element directory works unchanged.
- It pins `OPENBLAS/OMP/JULIA` threads internally and shells to
  `compute_bands.jl` (one SCF + bands + real-space `|ψ|²` binning).
- Determinism: the pinned DFT setup is deterministic, so the QP depths above
  reproduce to ~0.001 eV across runs.

---

## 6. What this probe shows about L3

- **L3 elicits real derivation, not recall:** 6 agents → 6 distinct closed forms,
  none the gold's `F_c` sine-transform.
- **The concealed test set does its job:** every derivation looked gold-quality on
  the visible Na/Al; the hidden K (diffuse-core alkali, strongest correction)
  exposed all of them. r5's K is *off by ~2×* despite a perfect-looking Na.
- **The tight per-element bars do their job:** the cruder-but-plausible models
  (which a loose 0.30 bar would have passed) are correctly held out.
- **The hard physics is the state-/momentum-resolved form factor** `f_c(|k+G|)` —
  the one ingredient none of the six reproduced, and the one needed to get K and Mg
  right simultaneously.
