# L3 probe — agent "ex2": the first passing solution

`ex2` is the **first agent (of 17) to pass Level 3**. It was run on the rewritten
L3 prompt (the paper's Eq (1)/(2) given; only `u_c` in the data) **after** two
requirements were added to `SETUP_L3` — *do a real step-by-step integrate-out (don't
posit/guess the coupling)* and *be careful with exchange*. Those two requirements are
what produced the passing route.

Solution preserved here: [`solution/run_qp.py`](solution/run_qp.py),
[`solution/qp_physics.py`](solution/qp_physics.py) (the closed form),
[`solution/dft_bands.jl`](solution/dft_bands.jl) (one pinned SCF + bands),
[`solution/method.md`](solution/method.md) (ex2's own derivation, verbatim).

---

## TL;DR

- **PASS.** Hidden K **0.135** (bar 0.17) and Mg **0.193** (bar 0.21) — both clear,
  all points covered, well below the KS baseline (K 0.614, Mg 0.434), errors
  unbiased. Not a flooding/shape artifact.
- It is the **first agent to derive the gold's *state-dependent momentum form
  factor*** — `z_core` built from `Σ_G c_nk(G)·(sine-transform of a core potential)`,
  not a scalar. The previous 16 all collapsed the coupling to a scalar and stayed
  single-`z`.
- It is **not** the exact gold: it weights the form factor by the **full core
  Hartree potential `V_H^core`** instead of the **orthogonalised per-channel
  fluctuation `(V_H_c − J_c)`** — i.e. it is missing the `−J_c` subtraction (plus a
  combinatorial factor 2). The structure is right; the potential is an approximation.
- The decisive ingredient is **state-dependence**, not the exact potential: even with
  the wrong (un-orthogonalised) weight, the `Σ_G c_nk(G)` contraction makes K drop
  from ~0.32 (single-`z`) to 0.135.

## Results

### Hidden metals (the benchmark verdict)

| element | z_core | KS RMSE | **QP RMSE** | bar | verdict |
|---------|--------|--------:|------------:|-----|---------|
| **K** (hidden)  | state-dep | 0.614 | **0.135** | 0.17 | PASS |
| **Mg** (hidden) | state-dep | 0.434 | **0.193** | 0.21 | PASS |
| **overall** | | | **mean 0.164** | | **PASS** |

K: `n_scored=10`, Mg: `n_scored=88`; both `n_missing=0`, `n_over_band_cap=0`.

### Public dev (Na, Al) and band-bottom depths

| element | z_core (band bottom) | KS RMSE | QP RMSE | Γ-depth KS→QP | reference |
|---------|----------------------|--------:|--------:|----------------|-----------|
| Na | ≈0.89 | 0.41 | **0.18** | −3.26 → −2.90 | expt −2.65…−2.78; eDMFT −2.84 |
| Al | ≈0.94 | 0.41 | **0.23** | −11.18 → −10.55 | ARPES bottom ≈ −10.58 |

Note the **dev-set irony**: ex2's full-band Na RMSE (0.18) is *worse* than the
single-`z` agents (~0.08). On the public data alone it would rank near the bottom —
yet it is the only one that **passes the hidden set**. State-dependence beats a
dev-tuned single-`z`; the dev set alone is a misleading ranker (cf.
`L3_prompt_probes.md`).

## The derivation (faithful summary)

Full text in [`solution/method.md`](solution/method.md). The closed form:

```
z_core_{nu,k} = 1 / ( 1 + Σ_c |g_{c,nk}|² / ΔE_c² )
g_{c,nk}      = 2 √(4π) ∫ u_c(r) · V_H^core(r) · ψ0(r) · r dr
ψ0(r)         = (1/√Ω_at) Σ_G c_nk(G) j₀(|k+G| r)          (Bloch s-wave at the atom)
V_H^core(r)   = monopole Hartree potential of the full core density (built from u_c)
ε_QP − ε_F    = z_core · (ε_KS − ε_F)
```

1. **Split the action; write the antisymmetrised vertex.** `ψ = ψ_v + ψ_c`; the
   valence–core Coulomb coupling has a **direct (Hartree)** piece `ρ_v ρ_c` and a
   **Fock (exchange)** piece (same fermion species). *(This is the step new1/new2/ex1
   skipped — ex2 actually wrote the interaction.)*
2. **Integrate the core out; split static/dynamical.** `g_v⁻¹ = g_0⁻¹ + δV_pp`;
   static part = `V_PSP` (already in `H_KS`), dynamical part = `z_core`.
3. **Closure + the even/odd-in-ω argument (correct, and decisive).** With the core
   collapsed to one mode at `ΔE_c`, `Σ_nk(ω) = Σ_c |g_{c,nk}|²/(ω − ΔE_c)`. A *direct*
   density–density response is **even in ω** → only a static shift, **no `iω` term**
   ("a direct picture misses `z_core` entirely"); the **exchange** resolvent
   `1/(ω−ΔE_c)` is **odd** → it gives the `iω` slope. Hence
   `z_core = [1 − ∂Σ/∂ω]⁻¹ = 1/(1 + Σ|g|²/ΔE²)`. *(new2 made the same even-in-ω
   observation but mis-used it; ex2 used it correctly.)*
4. **The coupling — state-dependent matrix element.** `g_{c,nk}` is the matrix element
   of the core potential between the core orbital `φ_c` and the Bloch state `ψ_nk`, so
   it enters through `ψ0 = Σ_G c_nk(G) j₀(|k+G|r)` → **state-dependent**.

### Why this is the gold's structure

Since `j₀(Kr) = sin(Kr)/(Kr)`, ex2's coupling is exactly
`Σ_G c_nk(G)·[√(4π)/K ∫ u_c · V_H^core · sin(Kr) dr]` — the gold's
`F_c(nk) = Σ_G c_nk(G) f_c(|k+G|)` with `f_c(K) = √(4π)/K ∫ u_c·(V_H_c − J_c)·sin(Kr)dr`.
**Same `Σ_G c_nk(G)·(sine-transform form factor)` structure.** The first of 17 agents
to reach it.

## The problem — where it falls short of the gold

| | gold | ex2 |
|---|---|---|
| weight in `f_c` | `(V_H_c − J_c)` (per-channel **fluctuation**, orthogonalised) | `V_H^core` (full core Hartree, **no `−J_c`**) |
| prefactor | exactly from the derivation | a **combinatorial factor 2** (closed-shell count) |

ex2 derived the **structure** (steps 1–3 are genuine: the antisymmetrised vertex, the
even/odd-in-ω split, the `z = 1/(1+Σ|g|²/ΔE²)` form, and the state-dependent matrix
element). But at step 4 it **did not evaluate the explicit 2nd-order Coulomb matrix
element**; it *substituted* a ready-made potential — "the core leaves behind its
Hartree field `V_H^core`, whose static part is `V_PSP`; its fluctuation is the
coupling." Computing the matrix element with valence–core antisymmetry/orthogonality
is what would have produced `(V_H_c − J_c)` (the `−J_c` is the orthogonalisation term);
substituting `V_H^core` skips it.

So ex2 is **half-derived**: structure derived from first principles, coupling
*operator* fixed by a physically-reasoned substitution rather than a computed
integral. (See `L3_prompt_probes.md` for why it stopped there: the explicit matrix
element is long; `V_H^core` is the obvious self-consistent object; and its band-bottom
matched eDMFT, so the dev self-test gave no signal that `−J_c` was missing — the `−J_c`
effect is invisible on the tight-core Na/Al and only matters on the diffuse-core K.)

## Why it passes anyway

The **state-dependent `Σ_G c_nk(G)` contraction** is the load-bearing ingredient.
The hidden K (diffuse 3s) has a band bottom with far more core overlap than its
near-`E_F` states; a single `z` cannot capture that (best single-`z` K ≈ 0.19–0.32),
but ex2's state-dependent form does — K = 0.135. The `−J_c` error shifts the
*magnitude* of the per-channel coupling but does not destroy the *state-dependence*,
so the dispersion shape is right enough to clear both bars. Headroom remains (it is
not the exact gold), which is healthy: the benchmark still separates a good
state-dependent approximation (passes) from the exact result.

## Prompt lineage

ex2 only became possible after `SETUP_L3` was given the paper's Eq (1)/(2) **and** the
two new requirements:
- *"do a real step-by-step integrate-out … do not posit a model self-energy or fix
  the coupling by dimensional analysis"* → forced the matrix element instead of a
  scalar (→ the form-factor structure);
- *"be careful with exchange … a purely direct density-response picture misses it"* →
  routed ex2 through the even-in-ω argument onto the exchange channel and the
  **state-dependent** Bloch-state matrix element.

Tally: **old / reframed / Eq-only prompts → 0 PASS / 14; this round → 1 PASS / 3.**

## Reproduction

`solution/run_qp.py` (CLI `--element-config --grid --out`) calls
`solution/dft_bands.jl` for one pinned LDA/GTH SCF + bands at `k = t·endpoint_frac`,
then applies `solution/qp_physics.py`'s `z_core_state` per Bloch state and writes
`E_pred_eV = z_core·(ε_KS − ε_F)`. Parameter-free, no per-element branch. The hidden
inputs ship only `u_c` + `ΔE_c`; `V_H^core` is built from `u_c` by the radial Poisson
solve in `qp_physics.radial_hartree`.
