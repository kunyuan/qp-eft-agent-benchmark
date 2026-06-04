# Benchmark Design: Frozen-Core Quasiparticle Bands

Status: **built & validated** (branch `benchmark-redesign`, PR #1).
Source physics: *Kohn–Sham Hamiltonian from Effective Field Theory: Quasiparticle
Band Narrowing from Frozen Core Dynamics* (arXiv:2604.25199).

This document is the design contract — what the benchmark measures, how it is
graded, and the decisions behind it. It is **not** agent-facing. As-built status:
three levels (L1/L2/L3) are assembled (`agent_packet/`, `evaluator/hidden/`) and
packaged as Harbor tasks (`harbor/`); the gold reproduces the paper; fresh solver
agents pass L1 (3 independent runs) and L2 through the hardened evaluator.

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
| **L2** | `z_core` closed form; atomic core radial data `u_c(r)`, `V_{H,c}(r)`, `DeltaE_c`; `J_c` computed from radial data | `f_c(K)` tables, `z_core_gamma` | implementing the form-factor quadrature (Eq. 7), extracting `c_nk(G)` from DFTK, assembling `z_core(n,k)` |
| **L3** | physical setup only (the two conditions, dual-fermion mechanism narrative, "missing physics is dynamical frozen-core"); atomic core radial data | the entire `z_core` / `F_c` formula | deriving `z_core` from the EFT and implementing it |

All three are scored identically on held-out elements; "generalization" is inherent
because the agent writes code against Na/Al only and is graded on concealed metals.

---

## 3. Evaluation protocol

- **Development set (agent-visible):** Na, Al.
- **Graded set (hidden):** K, Mg (Ca/Li/Si deferred — see §7 notes). Only
  `agent_packet/` (Na, Al) is given to the solver; the held-out metals are
  concealed.
- **DFTK in loop:** the evaluator runs the submitted `run_qp.py` on each hidden
  config+grid; the runner internally calls DFTK to produce KS bands, then applies
  the correction.
- **Ground truth:** the held-out ARPES references (this repo, `arpes_reference.csv`)
  are the falsifiable target. A **gold reference solution** (this repo, not
  agent-facing) reproduces the paper's Table I and supplies (a) the band
  cardinality at each k-point (occupied + first unoccupied) and (b) calibration of
  the RMSE thresholds.

### Anti-cheat (as built)

The element identity is **kept** in the config (the runner needs it to load the
pseudopotential, and `Z_nuclear` reveals it anyway — hiding it buys nothing). The
concealment is the test *set*, not the runtime name:

1. **Sanitized runner inputs** — the runner is handed an input dir with **no
   `arpes_reference.csv`**, so it cannot read the answers sitting next to the
   config (a real hole — the agent could discover it on Na/Al at dev time). In the
   Harbor tasks the runner also runs as `nobody` with `/tests/{hidden,gold}`
   root-only.
2. **Offline verifier** — the verifier phase runs with **no network**
   (`[verifier] network_mode = "no-network"` in each `task.toml`), so the submitted
   `run_qp.py` cannot fetch the held-out answers online (e.g. `git clone` the
   authors' `eft-psp` repo, which has the K/Mg band data). DFTK runs offline because
   the pinned env + GTH pseudo artifact are baked into the image at build. (The
   *agent* phase keeps network — it needs its own model API — but sees only Na/Al,
   and any network call its code makes fails at verify time. The override requires a
   Harbor provider that supports dynamic network policy.)
3. **Concealed test set + no-hardcode audit** — the solver writes generic code
   against Na/Al only; memorized per-element numbers can't enter code that must run
   on unknown metals without a forbidden per-element branch.
4. **KS-baseline gate (§4)** — the correction must be *necessary* to pass, so a
   memorized/plumbing-only submission fails.

---

## 4. Scoring (as built)

The original scorer picked, per `point_id`, the *nearest of an unlimited set* of
predicted bands → flooding the energy window drove RMSE → 0 with zero physics.
Fatal. The scorer (`evaluator/validate_submission.py`, mirrored in each Harbor
`tests/score.py`) now:

- **Band-set shape (`INVALID_SHAPE`):** every reference point must be predicted with
  **exactly** the gold band count per k-point (occupied bands + the first unoccupied
  band — see §5c). Closes flooding, "sparse" (drop hard points), and "under-band"
  cheats. Honest agents on the pinned setup match exactly.
- **Nearest-band RMSE** is then safe (n_pred = 1–4 bands, well-separated).
- **KS-baseline gate:** a bare uncorrected-KS submission scores ~0.4–0.6 eV and
  FAILs; `ks_baseline_rmse_eV` is reported for audit. The correction is *necessary*.
- **Sanitized inputs:** the runner never receives `arpes_reference.csv` (§3).
- **Thresholds** (calibrated from the gold): PASS < 0.30, PARTIAL 0.30–0.40, FAIL.
  Overall = aggregate over the hidden set; any per-element disqualifier
  (`REJECTED_FLOODING` / `INVALID_SHAPE` / `NO_PREDICTION`) sinks the submission.

---

## 5. Engineering (as built)

- **Structure spec:** the config carries the explicit `lattice_vectors_bohr`
  (columns of the DFTK lattice) and `atom_positions_frac` (≥1 atom). This is INPUT
  (like the pseudo), not the answer; it's required because the public dev set is
  cubic 1-atom (Na bcc, Al fcc) while hidden Mg is hcp 2-atom — without it a generic
  solver cannot build Mg's cell. Path: `endpoint_frac` + `k_frac = t*endpoint_frac`.
- **DFTK reproducibility:** DFTK `0.7.25`, pseudo family
  `cp2k.nc.sr.lda.v0_1.largecore.gth` (valence Na/K/Li q1, Mg/Ca q2, Al q3, Si q4),
  **per-element Ecut (15 for Na/K, 20 for Li/Ca/Mg/Al/Si)** and kgrid matching the
  authors' setup, Fermi–Dirac smearing `0.001 Ha`. Manifest committed under
  `environment/`.
- **Atomic core data schema (L2/L3):** per core s-channel, `atomic_core_<c>.csv`
  (`r_bohr,u_c,V_H_c`) + `core_model.json` (`DeltaE_c`). L1 instead ships
  `fc_table_<c>.csv` (precomputed `f_c(K)`). Generated by the gold's radial atomic
  solver; committed as data files.

## 5c. Which bands to emit (occupied + first unoccupied)

At each k-point, ascending in energy, emit every occupied band (`E_KS < E_F`) plus
the single lowest unoccupied band, then stop (count = `n_occ + 1`). No empirical eV
margin: this captures a Fermi-crossing band edge (where ARPES resolves a band just
below E_F that DFT places just above) and is robust to the exact εF placement at the
edge. (An earlier `+0.5 eV` margin was replaced because it was arbitrary and
εF-sensitive.) `E_pred = E_QP − E_F` may be slightly positive for that first
unoccupied band.

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

Notes (resolved during validation):
- **Multi-band output is mandatory** for Mg/Al (the ARPES tracks several dispersing
  bands). The runner emits occupied + first unoccupied (§5c).
- **Nearest-band matching is safe here**, not gameable: the flooding/`INVALID_SHAPE`
  guards force exactly the gold band count, and the occupied bands are 3–4 eV apart,
  so nearest-match is effectively per-band assignment. Verified on Mg: each ARPES
  point lands coherently on its band (deep→band1, mid→band2, shallow→band3) with
  small residuals — genuine per-band agreement, not loose matching. (A Hungarian
  assignment was considered and found unnecessary.)
- **Mg path is fine.** The earlier worry that codex's grids were mis-mapped was a
  misread (I'd taken a single deep ARPES point as Γ): the per-band check shows the
  Mg gold reproduces the experiment band-by-band, consistent with the paper. No grid
  regeneration was needed.
- **Al's worst point** is a single Γ-X Fermi-crossing band edge (DFT +0.08, ARPES
  −0.23 eV) — handled by the occupied+first-unoccupied rule (§5c). Al is public/dev,
  not scored.

## 6. Status & remaining

**Done:** gold reproduces Table I (Na KS −3.266 / QP −2.616); L1/L2/L3 packets +
hidden sets assembled; evaluator scoring (shape guard + KS gate + sanitized inputs)
+ Harbor tasks; band rule and structure spec finalized; fresh agents pass L1 (×3)
and L2 through the hardened evaluator; verifier hardening integrated (was codex
PR#2). Threshold/anti-cheat regression tests pass.

**Remaining / optional:**
- Container `run-as-nobody` needs a Docker run to confirm depot permissions (no
  Docker on the build host).
- Optional Li/Ca calibration anchors (Ca's EFT-vs-experiment gap is the worst, 0.24
  eV; Li has only eDMFT, no ARPES — frame as reference, not fit targets).
- Optional L3 agent run (frontier "derive the formula" rung; not yet agent-tested).

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
observable, different); LiH is a compound test. **As-built benchmark set: Na, Al
(public dev) + K, Mg (hidden).** Ca/Li/Si = optional future additions (Ca needs
EMS data digitized; Li has no ARPES; Si's correction is too small to discriminate).
