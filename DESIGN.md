# Benchmark Design: Frozen-Core Quasiparticle Bands

Status: **draft / in progress** (branch `benchmark-redesign`).
Source physics: *Kohn–Sham Hamiltonian from Effective Field Theory: Quasiparticle
Band Narrowing from Frozen Core Dynamics* (arXiv:2604.25199).

This document is the design contract. It records what the benchmark measures, how
it is graded, and the decisions behind it. It is **not** agent-facing.

---

## 1. What we are measuring

An agent's ability to do **first-principles theoretical-physics reasoning that
manifests as correct, generic code**. The single physics problem is the
frozen-core quasiparticle band-narrowing correction:

```
E_QP(n,k) - E_F = z_core(n,k) * (E_KS(n,k) - E_F)
z_core(n,k)     = 1 / (1 + sum_c |F_c(n,k)|^2 / DeltaE_c^2)
F_c(n,k)        = sum_G c_nk(G) * f_c(|k+G|)
```

The Kohn–Sham band structure is computed with **DFTK in the loop** (LDA / GTH /
plane waves). The agent must produce a generic `run_qp.py` that, given an element
config and a k-grid, returns occupied quasiparticle band energies.

### Core design principle (consequence of "only grade numbers")

We do **not** grade the derivation text. Theory reasoning is tested *implicitly*:
to get the numbers right the agent must have reconstructed the correct theory.
**This only works if the formula and numerical shortcuts are withheld** — otherwise
"numbers correct" ≠ "theory correct" and the task degrades to DFTK plumbing.
Withholding is therefore the mechanism that keeps reasoning in the loop.

---

## 2. Difficulty ladder (L1 / L2 / L3)

Same physics problem, same numerical scoring, same DFTK-in-loop eval protocol.
The levels differ only in **how much the `agent_packet/` withholds**.

| Level | Packet gives | Packet withholds | Primarily tests |
|-------|--------------|------------------|-----------------|
| **L1** | `z_core` closed form **and** precomputed `f_c(K)` tables + `DeltaE_c` | numerical answer (`z_core_gamma`), per-element values | DFTK wiring, correct application, generalization (floor) |
| **L2** | `z_core` closed form; atomic core radial data `u_c(r)`, `V_{H,c}(r)`, `J_c`, `DeltaE_c` | `f_c(K)` tables, `z_core_gamma` | implementing the form-factor quadrature (Eq. 7), extracting `c_nk(G)` from DFTK, assembling `z_core(n,k)` |
| **L3** | physical setup only (the two conditions, dual-fermion mechanism narrative, "missing physics is dynamical frozen-core"); atomic core radial data | the entire `z_core` / `F_c` formula | deriving `z_core` from the EFT and implementing it |

All three are scored identically on held-out elements; "generalization" is inherent
because the agent writes code against Na/Al only and is graded on concealed metals.

---

## 3. Evaluation protocol

- **Development set (agent-visible):** Na, Al.
- **Graded set (hidden, anonymized):** K, Mg, Ca, Li, Si (+ spares). The element
  identity is **not** revealed to the solver — only `agent_packet/` (Na, Al) is.
- **DFTK in loop:** the evaluator runs the submitted `run_qp.py` in an **offline
  sandbox** (no network) on each hidden config+grid; the runner internally calls
  DFTK to produce KS bands, then applies the correction.
- **Ground truth:** the held-out ARPES references (this repo, `arpes_reference.csv`)
  are the falsifiable target. A **gold reference solution** (this repo, not
  agent-facing) reproduces the paper's Table I and supplies (a) the occupied-band
  cardinality at each k-point and (b) calibration of the RMSE thresholds.

### Anti-cheat (decided: rely on no-network + concealed set, no counterfactuals)

The solver never sees the test elements when writing code, so memorized ARPES
numbers cannot be injected into code that must also work on unknown metals — unless
it hardcodes per-element branches, which are forbidden and cannot survive the
concealed set. Hardening:

1. Strip the `element` name field from hidden configs (the code must not need it).
2. No-network sandbox for the runner at eval time.
3. No-hardcode audit + the KS-baseline gate (§4): the correction must be *necessary*
   to pass, so a memorized/plumbing-only submission fails.

---

## 4. Scoring (fixes the current exploit)

**Current bug:** `validate_submission.py` picks, per `point_id`, the *nearest of an
unlimited set* of predicted bands. Flooding predictions across the energy window
drives RMSE → 0 with zero physics and no DFTK. Fatal.

**Fix:**
- The runner must output **exactly the occupied band set** at each k-point. The
  evaluator knows the occupied-band cardinality `n_occ(point_id)` from the gold
  reference; submissions with `n_pred != n_occ` are penalized (extra/missing bands).
- Match predicted ↔ reference bands by **one-to-one assignment** (Hungarian), not
  nearest-of-many.
- **KS-baseline gate:** a bare uncorrected-KS submission must score FAIL. Passing
  requires the correction to close the KS→ARPES gap. This is the real signal that
  the agent did the physics.
- Per-element verdict by RMSE thresholds (calibrated from the gold reference, not
  the current hand-set 0.20/0.40); overall = aggregate over the hidden set.

---

## 5. Engineering fixes (no decision needed)

- **k-path spec:** define `Gamma-N`, `Gamma-X`, `Gamma-A` as exact reciprocal
  fractional coordinates (Setyawan–Curtarolo) so two correct implementations agree
  on the k-vectors. Map grid `t=0 -> Gamma`, `t=1 -> endpoint`.
- **DFTK reproducibility:** pin DFTK `0.7.25`, pseudo family
  `cp2k.nc.sr.lda.v0_1.largecore.gth` (valence counts match the configs:
  Na/K/Li q1, Mg/Ca q2, Al q3, Si q4), converged `Ecut` (validate 30 Ha is enough),
  Fermi–Dirac smearing `0.001 Ha`. Manifest committed under `environment/`.
- **Atomic core data schema (L2/L3):** per element, the dominant core channel(s)
  `c` with reduced radial wavefunction `u_c(r)` on a uniform grid, orbital Hartree
  potential `V_{H,c}(r)`, self-Coulomb `J_c`, and excitation energy `DeltaE_c`.
  Generated by the gold reference's radial atomic solver; committed as data files.

---

## 5b. Calibration results (gold_runner.jl vs real ARPES) — VALIDATED

All four simple metals confirm the correction is **necessary and sufficient**:

| El | KS RMSE (eV) | QP RMSE (eV) | improvement | note |
|----|--------------|--------------|-------------|------|
| Na | 0.413 (FAIL) | 0.078 (PASS) | 5.3× | alkali, large correction |
| K  | 0.614 (FAIL) | 0.139 (PASS) | 4.4× | alkali, largest correction |
| Mg | 0.434 (FAIL) | 0.187 (PASS) | 2.3× | needs multi-band output |
| Al | 0.414 (FAIL) | 0.248 (PASS) | 1.7× | small correction, tightest margin |

**Calibrated thresholds: PASS < 0.30 eV, PARTIAL 0.30–0.40, FAIL > 0.40 eV.**
Separates every QP from every KS; KS-baseline gate satisfied automatically.

**Per-band before/after** (`reference/validate_bands.py`, same band assignment):
the correction acts exactly where the physics says it should. For multi-band
Mg/Al the *deep* band carries the DFT overbinding (KS RMSE 0.42–0.64 eV, mean
−0.4 to −0.6) and is pulled into agreement (2–3×), while bands *near E_F* barely
move (1.0×) because `z_core·(E_KS−E_F) → 0` at the Fermi level. The systematic
overbinding (negative KS mean) is removed in every case (QP mean → ~0). The
coherent per-band match (deep→band1, mid→band2, shallow→band3 for Mg, with the
experiment tracking the real bands) also confirms the k-mapping is sound for the
multi-band elements — the score is genuine per-band agreement, not loose
nearest-of-many matching.

| El | DFT band-1 RMSE / mean | QP band-1 RMSE / mean |
|----|------------------------|-----------------------|
| Na | 0.413 / −0.366 | 0.079 / −0.009 |
| K  | 0.614 / −0.599 | 0.139 / −0.033 |
| Mg | 0.635 / −0.587 | 0.213 / −0.029 |
| Al | 0.423 / −0.390 | 0.231 / −0.130 |

Caveats found during validation (feed into the evaluator + packet rebuild):
- **Multi-band output is mandatory** for Mg/Al/Ca/Si (the ARPES tracks a band that
  disperses; lowest-band-only fails). The runner outputs all occupied bands.
- **Nearest-band matching is forgiving** (it's how validation passed) but is the
  gameable path — the real evaluator must use one-to-one assignment vs the gold
  occupied-band set, not nearest-of-many.
- **Per-element path mapping is inconsistent in codex's grids.** Each experiment
  uses different x-units (Na: Å⁻¹ ÷ |Γ-N|; K: already fractional; Mg: Γ-A-Γ folded,
  x extends to 1.67×|Γ-A|; Al: Γ-X). The benchmark grids + ARPES references must be
  **regenerated from the raw `reference/expts/*.csv`** with one correct convention,
  not trusted from codex.

## 6. Open build tasks

1. **[in progress]** Validate DFTK KS bandwidth for Na (`E(Γ)-E_F ≈ -3.27 eV`).
2. Gold reference `run_qp.py`: DFTK valence bands + atomic `f_c` solver +
   `z_core(n,k)` → reproduce Table I for Na/Al/K/Mg/Ca/Li/Si. Calibrate thresholds.
3. Fix evaluator scoring + KS-baseline gate using gold cardinalities.
4. Restructure `agent_packet/` into L1/L2/L3; remove `z_core_gamma` and `element`
   from hidden configs; add k-path spec and atomic core data files.

---

## 7. Pinned facts (for reproducibility) — RESOLVED against the authors' repo

Ground truth is the authors' pipeline `github.com/iintSjds/eft-psp` (two stages:
`<el>/run_ks.jl` → DFTK KS data; `<el>/freq_correction.jl` + `atomic_hf.jl` →
atomic solver + coherent z_core correction). **Reproduced exactly in our env**
(Na: KS Γ −3.266 eV, QP Γ −2.616 eV, ratio 0.801).

- DFTK 0.7.25, PseudoPotentialData 0.3.2, Julia 1.12.1.
- Pseudo family: `cp2k.nc.sr.lda.v0_1.largecore.gth` (Zion: Na/K/Li 1, Mg/Ca 2, Al 3, Si 4).
- **Lattice constants are EXPERIMENTAL, not LDA-relaxed** (the earlier "LDA-relaxed"
  plan was wrong — the large-core PSP's LDA equilibrium is anomalous and is NOT what
  the paper used). The 3.27-vs-3.145(free-electron) gap is set by the paper's **low
  per-element Ecut**, not the lattice constant.
- Atomic solver: radial LDA with **Dirac exchange only** (no correlation), uniform
  grid dr=0.002, r_max=40 Bohr. `ΔE_c` from **ΔSCF** total-energy differences
  (hole in core orbital c), summed over all core s-channels.
- Form factor (paper Eq. 7, s-wave core):
  `f_c^K = (sqrt(4π)/K) ∫ u_c(r)[V_{H,c}(r) - J_c] sin(Kr) dr`, single-orbital
  Hartree `V_{H,c}`, `J_c = ∫ u_c^2 V_{H,c} dr`.
- Correction: `Δ(n,k) = Σ_c |Σ_G c_nk(G) f_c(|k+G|)/ΔE_c|^2`,
  `E_QP = εF + (E_KS - εF)/(1+Δ)`  (i.e. `z_core = 1/(1+Δ)`).

### Per-element gold settings + ground-truth (from the authors' `*_summary.txt`)

| El | struct | a (Bohr) | Ecut | kgrid | Zval | dom ch | ΔE_c (Ha) | KS Γ (eV) | QP Γ (eV) | ratio |
|----|--------|----------|------|-------|------|--------|-----------|-----------|-----------|-------|
| Li | bcc    | 6.632    | 20   | 8³    | 1    | 1s     | 2.7715    | −3.483    | −2.616    | 0.751 |
| Na | bcc    | 8.107    | 15   | 8³    | 1    | 2s     | 2.7006    | −3.266    | −2.616    | 0.801 |
| K  | bcc    | 9.874    | 15   | 8³    | 1    | 3s     | 1.7184    | −2.270    | −1.492    | 0.657 |
| Ca | fcc    | 10.545   | 20   | 8³    | 2    | 3s     | 2.5071    | −3.700    | −3.057    | 0.826 |
| Mg | hcp    | 6.066 (c/a 1.624) | 20 | 8×8×6 | 2 | 2s | 4.0662 | −6.921 | −6.310 | 0.912 |
| Al | fcc    | 7.653    | 20   | 8³    | 3    | 2s     | 5.6901    | −11.179   | −10.704   | 0.957 |
| Si | diamond | 10.263  | 20   | 6³    | 4    | 2s     | —         | (gap)     | —         | —     |

Notes: Li uses an analytic hydrogenic 1s form factor; Si is a semiconductor (gap
observable, different); LiH is a compound test. **Core simple-metal set for the
benchmark: Na, Al (public dev) + K, Mg, Ca, Li (hidden).** Si/LiH = advanced/optional.
