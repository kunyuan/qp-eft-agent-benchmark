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
| **L4** | the problem (puzzle + all-electron action + double-counting fact), the public ARPES, **Li as a closed-form development element** (config+grid, theoretical LDA/eDMFT anchor, mandatory Stage 1: complete analytic derivation + contraction enumeration on the two-electron core before generalizing); gold's public-element scores stated as the baseline | everything else: formula, structural ansatz, ALL atomic data, prescribed approximations (closure / s-channel are no longer commanded) | open frontier — own atomic solver, every truncation declared+controlled (consistency ledger in method.md), and *beyond-leading-order* physics: `beats_gold` / `beats_gold_all` reported by the evaluator. PASS bars unchanged (= leading-order quality); beating gold on every hidden element is the scientific result |

All levels are scored identically on held-out elements; "generalization" is inherent
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
2. **Concealed test set + no-hardcode audit** — the solver writes generic code
   against Na/Al only; memorized per-element numbers can't enter code that must run
   on unknown metals without a forbidden per-element branch.
3. **No network / no external lookup** — two instruction rules: (a) `run_qp.py`
   must not reach the network at run time (compute from local inputs + DFTK,
   offline); (b) the agent must not consult **any** external source while solving —
   no web search of any kind: not the underlying paper / derivation / answers / the
   `eft-psp` repo, and **not even library/API documentation**. (b) matters most for
   **L3** (the public paper would hand over the derivation), but applies throughout.
   This is solvable without web because DFTK is installed locally — the intended way
   to learn its API (e.g. `ExplicitKpoints`, `Gplusk_vectors_cart`) is Julia
   introspection of the installed package (`names`, `?`, `methods`, source), not the
   web. (Reflects how the reference solution was found.)
   - **Enforcement (rule → control):** these are rules by default. To enforce:
     run the verifier offline (Harbor `[verifier] network_mode = "no-network"`, or
     a `--network none` sandbox), and restrict the *agent* phase to an allowlist of
     only its model-API host (Harbor `[agent] network_mode = "allowlist"`) or run it
     with no web tools. Both are provider/harness-dependent and **not** set by
     default (a per-phase network override needs a provider with dynamic policy,
     e.g. E2B).
   - **Residual limit (cannot be controlled):** the paper is public and likely in
     the agent's *training data*, so a strong model could recall the closed-form
     `z_core` even fully offline. Network controls stop *live* lookup, not memory.
     This is intrinsic to **L3** (derive a published result); L1/L2 are unaffected
     (the formula is given there, and the hidden answers are per-point band data
     that can't be memorized precisely and aren't known to be the test set).
4. **Vertex diagonal audit (L3, post-hoc)** — for the submission's derived
   per-channel coupling potential, `∫ u_c² V_vertex dr` must be ~0 (the dynamical
   vertex is a zero-mean fluctuation against the core orbital; its mean is already
   inside `DeltaE_c` and the static PSP). A large nonzero diagonal flags a coupling
   selected by numerical coincidence on Na/Al rather than derived — the dominant
   failure mode of the strongest probed agents (see `evaluator/README.md` checklist,
   `experiments/L3_fable5/REPORT.md`, `maintainer_sources/NOTE-vertex-derivation-Li.md`).
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
- **Per-element thresholds (maintainer-only), overall = all-must-pass.** PASS for an
  element iff `rmse < bar(element)`; PARTIAL within `+0.10`; else FAIL. Overall PASS
  iff *every* hidden element clears its own bar (no averaging away a miss); any
  disqualifier (`REJECTED_FLOODING`/`INVALID_SHAPE`/`NO_PREDICTION`) sinks it. The
  bars are calibrated tight to the gold's hidden-set RMSE + a small margin:
  `K < 0.17`, `Mg < 0.21` (gold 0.139 / 0.187; faithful implementations reproduce
  the gold to ~0.001 eV across independent runs incl. no-web, so the margin is
  ~10–20× the observed spread). Why tight: the old global 0.30 was pinned by the
  gold's worst case **Al (0.248)** — but Al is a *dev* element, not hidden, so the
  hidden K/Mg bar can sit right at their gold. This correctly fails a cruder
  approximation that gets the magnitude but mis-fits the band k-dependence (e.g. a
  state-independent single-Z model, ~0.188 / 0.220 → PARTIAL, not PASS), which a
  loose 0.30 let through.

### 4b. Accuracy target stated to the agent (no grading number revealed)

The instruction does **not** state the pass bar. Instead it gives a *physical*
target: bare KS overshoots the bandwidth by 20–35%; a correct frozen-core
correction reaches ~0.1 eV agreement with experiment — the level of the many-body
reference eDMFT. Public calibration anchors are given as Γ-point depths (eV):

| el | LDA (Mandal '22) | eDMFT (Mandal '22) | expt |
|----|------|-------|------|
| Li | 3.48 | 2.60 | — |
| Na | 3.30 | 2.84 | 2.65–2.78 |
| Ca | 3.98 | 3.24 | 3.30 |

Only **Li/Na/Ca** appear (the hidden K/Mg are never shown). Caveats baked into the
wording: (i) the LDA column is literature — the agent uses its *own* pinned-DFTK KS
(matches to ~0.1 eV), it must not tune to these; (ii) "match eDMFT" is framed as
"~0.1 eV, comparable to eDMFT/experiment" — eDMFT itself is ~0.1 eV off experiment
(for K it is *further* from ARPES than the EFT method), so it is a peer reference,
not a higher-precision ground truth. The Γ-depth anchors pin the per-element
*magnitude* but not the k-dependence — so they help the agent calibrate, while the
real discriminator stays the hidden full-band RMSE with the tight per-element bars.

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

**Grading thresholds: tight per-element bars (hidden K < 0.17, Mg < 0.21 eV; §4).**
These sit just above the gold's hidden-set RMSE and well below every KS baseline, so
the KS-baseline gate is satisfied automatically and a cruder approximation that only
roughly matches the magnitude does not pass. (The public Al case, gold 0.248, is why
a *global* 0.30 was loose — but Al is dev, not hidden.)

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
