# L3 probe — agent "clean2": a single-`z` derivation, and what it reveals

`clean2` is one of two fresh, no-web agents run on the **reframed** L3 setup (the
EFT "integrate out the core" framing; see PR for `reference/build_packet.py`
`SETUP_L3`). It produced a clean, self-consistent derivation of a **single
quasiparticle weight `z` per element** and turned out to be — on the *bandwidth*
metric the literature reports — **as accurate as eDMFT and better than the
published EFT**, while still failing the benchmark's full-band test. That gap is
the point of this note.

Solution preserved here: [`solution/run_qp.py`](solution/run_qp.py),
[`solution/_ks_bands.jl`](solution/_ks_bands.jl),
[`solution/method.md`](solution/method.md) (clean2's own derivation, verbatim).
Predictions: [`outputs/`](outputs/). Ca inputs used for the cross-material
comparison (generated with the gold atomic solver): [`ca_case/`](ca_case/).

---

## TL;DR

- **Derivation (clean2):** integrate the core out → the static reduction is the
  pseudopotential DFTK uses; the leading missing piece is a 2nd-order, closure
  self-energy whose *slope* gives a quasiparticle weight
  **`z = 1 − λ`, `λ = Σ_c ⟨u_c² V_H_c⟩ / ΔE_c²`**, and `E_QP − E_F = z·(E_KS − E_F)`.
- **It is a single `z` per element (k-independent).** clean2 argues this from
  "high-energy on-site core fluctuation → local self-energy → the eDMFT picture."
- **`J_c` is not a controlled reduction of the gold's `f_c`** — its coupling is
  ~1.6× too weak; the accuracy comes from a *coincidental cancellation* against the
  `z = 1−λ` linearization (§2), with the largest residual on K (where it fails).
- **On the Γ-depth (occupied bandwidth) it is excellent** — mean |deviation from
  ARPES| **0.09 eV** across Na/Al/K/Ca/Mg, vs **0.14** for the published EFT (gold)
  and **0.10** for eDMFT. It is the *best of the three* on bandwidth.
- **But it fails the benchmark's full-band RMSE** (hidden K 0.190, Mg 0.211; bars
  0.17/0.21 → PARTIAL, not pass). The single `z` reproduces the band *bottom* but
  not the *dispersion* — the k-dependence only the state-resolved gold captures.
- **Lesson:** a single-`z` model can match (even beat) the published *bandwidth*
  numbers while being physically incomplete. The benchmark scores the **full ARPES
  band**, so it is strictly more discriminating than the one-number bandwidth
  comparison in the literature.

---

## 1. The derivation (faithful summary)

Full text in [`solution/method.md`](solution/method.md).

1. **Integrate out the core.** The all-electron problem reduces to an effective
   valence theory; its static part is the conventional pseudopotential (= the
   `H_KS` DFTK produces). The missing piece is the core's frequency-dependent
   response (a valence electron virtually excites a core particle–hole pair at
   `ΔE_c`) — the leading 2nd-order core self-energy.
2. **Closure** at the single high scale `ΔE_c` collapses each core s-channel to one
   bosonic mode `Σ_c(ω) = W_c/(ω − E_F + ΔE_c)`, whose integrated weight is the core
   orbital's Hartree self-energy `W_c = ⟨u_c² V_H_c⟩ = ∫ u_c(r)² V_H_c(r) dr`.
3. **Slope, not value.** Expanding about `E_F`: the constant `W_c/ΔE_c` is the
   static part already in the pseudopotential (cancels in `E − E_F`); the slope
   `Σ_c'(E_F) = −W_c/ΔE_c²` is the dynamical content, giving the quasiparticle
   weight

   ```
   z = 1 − λ,   λ = Σ_c ⟨u_c² V_H_c⟩ / ΔE_c²,   E_QP − E_F = z·(E_KS − E_F).
   ```

4. **Single `z` per element (clean2's claim).** "On-site, high-energy core
   fluctuation → local self-energy → k-independent weight — the eDMFT picture." The
   `1/ΔE_c²` denominator makes the deep 1s negligible and lets the soft outer core
   dominate, so the element trend rides mainly on `ΔE_{2s/3s}`.

Parameter-free: `⟨u_c² V_H_c⟩` and `ΔE_c` are atomic inputs; no fit, no per-element
branch. (This is the `J_c/ΔE_c²` family; cf. old probe r6.)

> **Where it is incomplete vs the gold.** The gold keeps the **state-resolved form
> factor** `Δ(nk)=Σ_c|Σ_G c_nk(G) f_c(|k+G|)|²/ΔE_c²` with
> `f_c(K)=√(4π)/K∫u_c(V_H_c−J_c)sin(Kr)dr`. clean2 drops the `c_nk(G)` overlap and
> the sine-transform form factor, replacing the state-resolved coupling by the
> single scalar `⟨u_c²V_H_c⟩` — i.e. it assumes the renormalization is k-independent.

---

## 2. Why `J_c` still lands accurate — a coincidental cancellation

A natural worry: clean2's coupling `J_c = ∫u_c² V_H_c` is **not** the gold's
coupling. In the form factor `f_c(K)=√(4π)/K∫u_c(V_H_c−J_c)sin(Kr)dr`, the coupling
is the **fluctuation** `(V_H_c − J_c)`; `J_c` is the constant **subtracted** to make
that fluctuation zero-mean (`∫u_c²(V_H_c−J_c)=0`). clean2 uses the *subtracted
constant itself* as the coupling. So how does it land within ~0.1 eV of the gold?
Two errors that partially cancel.

**(i) The coupling really is weaker — by ~1.6×, not 2×.** Compare `J_c` to the
gold's *band-bottom* effective coupling `|F_c(Γ)|² = (1/z_gold − 1)·ΔE_c²` — **not**
to `f_c(0)²` (the K→0 maximum, which over-states it, because the band-bottom Bloch
state is not a pure plane wave, `c_Γ(0)² ≈ 0.86`):

| el (dom ch) | gold band-bottom coupling | `f_c(0)²` | **`J_c`** | `J_c` / gold |
|-------------|--------------------------:|----------:|----------:|-------------:|
| Na 2s | 1.81 | 2.10 | 1.15 | 0.63 |
| K 3s  | 1.54 | 1.71 | 0.70 | 0.46 |
| Ca 3s | 1.32 | 1.41 | 0.77 | 0.58 |
| Mg 2s | 1.60 | 1.68 | 1.29 | 0.81 |

So `J_c` is **46–81 %** of the real coupling (mean ~0.6). (My first pass wrongly
compared `J_c` to `f_c(0)²` and quoted "~2×"; the band-bottom value is the right
reference.)

**(ii) A second error over-corrects and cancels most of it.** clean2's `z = 1 − λ`
is a *linearization* of `z = 1/(1+λ)`. For any `λ > 0`, `1 − λ < 1/(1+λ)`, so the
linearization **over-corrects** (narrows more than the true weight for the same
coupling). The weaker coupling pushes `z` up (under-correct); the linearization
pushes it down (over-correct); they partly cancel:

| el | gold z | clean2 `z = 1−λ` | if it used `1/(1+λ)` | residual |
|----|-------:|-----------------:|---------------------:|---------:|
| Na | 0.801 | 0.838 | 0.861 | +0.037 |
| K  | 0.657 | 0.748 | 0.799 | +0.091 |
| Ca | 0.826 | 0.868 | 0.883 | +0.042 |
| Mg | 0.912 | 0.919 | 0.925 | +0.007 |

Read the "`1/(1+λ)`" column as the weak coupling showing through (Na would be 0.861,
well above gold 0.801 — under-corrected); the `1 − λ` linearization pulls it back.
The two roughly cancel, leaving a small **net under-correction** (+0.007 … +0.091).

**The residual is largest for K** (+0.091): K's `J_c` is the most deficient (46 %)
*and* its `λ` is the largest (0.25), so the linearization helps but cannot fully
close the gap — which is exactly why K is clean2's worst element (RMSE 0.190, the
one furthest over its bar).

**So clean2's accuracy is not a controlled reduction of the gold.** `J_c` is not
`f_c`; the agreement comes from a fortunate cancellation of a ~1.6× weak coupling
against a linearization that over-corrects by a similar amount, with the `1/ΔE_c²`
factor carrying the dominant element trend. A *correct* single-`z` reduction would
evaluate the real `(V_H_c − J_c)` form factor at the band-bottom momenta — one
number per element, but the right one — not substitute `J_c`.

---

## 3. Results

### Full-band RMSE vs ARPES (the benchmark metric)

| element | λ | z | KS RMSE | **QP RMSE** | verdict |
|---------|--:|--:|--------:|------------:|---------|
| Na (dev)    | 0.162 | 0.838 | 0.413 | **0.089** | — |
| Al (dev)    | 0.047 | 0.953 | 0.414 | **0.209** | — |
| **K** (hidden)  | — | — | 0.614 | **0.190** | PARTIAL (bar 0.17) |
| **Mg** (hidden) | — | — | 0.434 | **0.211** | PARTIAL (bar 0.21) |
| **overall** | | | | | **PARTIAL — not pass** |

K is 0.02 over its bar; Mg is 0.001 over — the closest any L3 agent has come, but
still failing because a single `z` cannot reproduce the band dispersion.

### Cross-material bandwidth (Γ-point occupied depth, eV)

Ca was generated for this comparison with the gold atomic solver (gold reproduces
the eft-psp value 3.057 exactly, confirming the setup). All values are the band
bottom `−min(E_pred)`:

| el | LDA | eDMFT | **gold (EFT)** | **clean2 (1-`z`)** | ARPES |
|----|----:|------:|---------------:|-------------------:|------:|
| Na | 3.27 | 2.84 | 2.616 | **2.732** | 2.65–2.78 |
| Al | 11.18 | — | 10.704 | **10.656** | ~10.6 |
| K  | 2.27 | 1.42 | 1.492 | **1.693** | 1.60 |
| Ca | 3.70 | 3.24 | 3.057 | **3.209** | 3.30 |
| Mg | 6.92 | 6.18 | 6.310 | **6.359** | 6.15 |

### |deviation from ARPES| — lower is better

| el | eDMFT | gold (EFT) | **clean2 (1-`z`)** |
|----|------:|-----------:|-------------------:|
| Na | 0.12 | 0.10 | **0.02** |
| Al | — | 0.10 | **0.06** |
| K  | 0.18 | 0.11 | **0.09** |
| Ca | 0.06 | 0.24 | **0.09** |
| Mg | 0.03 | 0.16 | 0.21 |
| **mean** | **0.10** | **0.14** | **0.09** |

On the band-bottom depth, clean2's single `z` is **the most accurate of the three**
(it is closer to ARPES than the published EFT on K and Ca, where the EFT
under-corrects). Mg is its one weak point.

---

## 4. The finding: bandwidth ≠ band

clean2 is a near-perfect illustration of why the benchmark scores the **full ARPES
band**, not a single bandwidth number:

- **On the Γ-depth** (what the paper / eDMFT / G₀W₀ tables report), clean2's single
  `z` is competitive with eDMFT and beats the published EFT — a single scalar per
  element, no state-dependence, looks "publication-grade."
- **On the full band**, the single `z` rescales every k-point by the same factor,
  so it nails the band bottom (Γ, where the deviation is largest and most weighted
  in a one-number comparison) but mis-fits the intermediate dispersion. The
  benchmark's per-point RMSE over the whole occupied band catches exactly this, and
  the tight bars (K<0.17, Mg<0.21 — the latter sits right at the single-`z` ceiling)
  correctly hold it out.

So the benchmark is **strictly more discriminating than the literature bandwidth
comparison**: a model that would "look as good as eDMFT" in Table-I style still
fails here unless it carries the state-resolved form factor.

---

## 5. Reproduction

From the repo root, with Julia + the pinned project (`environment/`):

```bash
export JULIA_PROJECT="$PWD/environment"
export OPENBLAS_NUM_THREADS=4

# (a) run clean2 on a benchmark element (public Na/Al, hidden K/Mg)
python3 experiments/L3_clean2/solution/run_qp.py \
  --element-config evaluator/hidden/L3/K/element_config.json \
  --grid          evaluator/hidden/L3/K/grid.csv \
  --out           /tmp/clean2_K.csv
# QP Γ-depth = -min(E_pred_eV); full-band RMSE via evaluator/validate_submission.py --level 3

# (b) run clean2 on Ca (inputs vendored here; generated with the gold atomic solver)
python3 experiments/L3_clean2/solution/run_qp.py \
  --element-config experiments/L3_clean2/ca_case/element_config.json \
  --grid          experiments/L3_clean2/ca_case/grid.csv \
  --out           /tmp/clean2_Ca.csv
```

`run_qp.py` reads `core_model.json` + `atomic_core_*.csv` relative to the config,
runs one pinned SCF (`_ks_bands.jl`), and emits occupied + first-unoccupied bands
as `E_pred_eV = z·(E_KS − E_F)`. Expected QP Γ-depths: Na 2.732, Al 10.656,
K 1.693, Ca 3.209, Mg 6.359 eV (reproduce [`outputs/`](outputs/)).

Not covered: **Li** (eft-psp uses a special analytic 1s form factor) and **Si**
(semiconductor) — neither is in the gold solver's element set.
