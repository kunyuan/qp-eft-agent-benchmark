# Reference solution & validation (maintainer-only — NOT agent-facing)

This directory holds the **gold reference** for the benchmark and the validation
artifacts proving it is solvable and correctly calibrated. None of it is given to
a problem-solving agent.

## Provenance

The authoritative physics pipeline comes from the authors' research repo
**`github.com/iintSjds/eft-psp`** (the code behind arXiv:2604.25199):

- `gold/sodium/{run_ks.jl, freq_correction.jl, atomic_hf.jl}` — vetted copies of
  the authors' two-stage pipeline (DFTK KS → atomic ΔSCF + coherent z_core).
- `expts/*.csv` — the real experimental ARPES band data (`naband1/2`, `alband`,
  `kband`, `mgband`) from that repo's `notes/expts/`. These are the falsifiable
  scoring targets.

## Our additions

- `gold_runner.jl` — parametrized reference: takes an element's pinned DFT setup +
  a grid of fractional path coordinates `t`, computes KS bands at `k(t)=t·endpoint`,
  applies the **coherent** correction `Δ=Σ_c|Σ_G c(G)f_c(|k+G|)/ΔE_c|²`,
  `E_QP=εF+(E_KS−εF)/(1+Δ)`, and writes predictions through the benchmark interface.
- `atomic_solver.py` — an independent from-scratch reimplementation of the atomic
  radial LDA (Dirac-exchange) solver + form factor, used to cross-validate the gold
  (matches to 4 digits: Na `f(0)=1.451`, `Δ=0.289`).
- `na_check.jl` — minimal DFTK Na band sanity check.

## Validation status (Na)

```
Gold reproduces paper Table I:   KS Γ −3.266 eV,  QP Γ −2.616 eV,  z=0.801
Gold vs real ARPES (31 points):  KS RMSE 0.413 eV (FAIL)  →  QP RMSE 0.078 eV (PASS)
```

The correction is both **necessary** (bare KS fails) and **sufficient** (QP passes),
which is the property the benchmark scores. See `../DESIGN.md` for the full design.

## Running

```bash
ENV=environment   # pinned Julia project (DFTK 0.7.25, PseudoPotentialData 0.3.2)
julia --project=$ENV reference/gold_runner.jl Na agent_packet/data/public/Na/grid.csv /tmp/na_pred.csv
```

Generated `*.jld2` / `*.jls` / `*_summary.txt` are git-ignored (regenerate by running).
