# Level 2 — compute the form factor

The formula is given, but you must compute the core form factor `f_c(K)` yourself from the atomic core data. Read `THEORY.md`.

Public development elements: `Na/`, `Al/` (each with `element_config.json`, `grid.csv`, `arpes_reference.csv`, and the level's data files).

## Environment & how to compute the Kohn-Sham band

The compute environment provides **Julia + DFTK 0.7.25 + PseudoPotentialData
0.3.2** (and Python). `run_qp.py` is a Python entry point; it may shell out to a
Julia/DFTK script (the natural approach) or drive the DFT however you like, as
long as the KS step uses the pinned setup below.

The pseudopotential is a PseudoPotentialData family identifier — load it by name:

```julia
using DFTK, PseudoPotentialData
psp = load_psp(PseudoFamily("cp2k.nc.sr.lda.v0_1.largecore.gth"), Symbol(element))
```

From the DFTK run you need, per k-point: the KS eigenvalues, the Fermi level
`εF`, the plane-wave coefficients `c_nk(G)` of each Bloch state, the integer
`G`-vectors, and the reciprocal lattice (to form `|k+G|` in Bohr^-1). Compute
the bands at the explicit k-points `k_frac = t * endpoint_frac` (not a
density-based auto path), so each grid row maps to one k-point.

## Pinned DFT setup (use exactly; do not "improve" it)

LDA / GTH `cp2k.nc.sr.lda.v0_1.largecore.gth` / `Ecut_Ha`, `kgrid`,
`smearing_Ha` from the config / Fermi-Dirac smearing. The Kohn-Sham band is the
starting point — a converged KS band alone is NOT a quasiparticle prediction.

## Submission

Submit a directory containing:

- `run_qp.py` with this exact interface:

```bash
python run_qp.py --element-config element_config.json --grid grid.csv --out qp_bands.csv
```

- `method.md`: a short note on how you derived (Level 3) / implemented the
  correction. **Not scored** — read only for the maintainer audit (parameter-free,
  no per-element hardcoding). For Level 3, show your derivation.

### Output (`qp_bands.csv`)

Columns `element,point_id,t,E_pred_eV`, energies **relative to the Fermi level**.
`t` is the fractional path coordinate: `k_frac = t * endpoint_frac` (both in
`element_config.json`).

Output the **occupied** quasiparticle bands at each grid point — one row per band
with `E_pred_eV < 0` (below E_F). Compute enough bands (a few above `z_valence`)
to cover all occupied states, then keep only the occupied ones; do **not** emit
unoccupied bands. The evaluator caps predictions to the true occupied-band count
per point, so extra/flooded bands are rejected.

## Rules

- Parameter-free: no fitting to ARPES, no per-element tuning, no hardcoded
  output values, no per-element `if` branches. The SAME code path runs on the
  public elements and on concealed held-out metals.
- Develop only against the public elements here (Na, Al). Hidden metals are
  scored by the evaluator and are not revealed.

## Scoring

Predictions are compared to held-out ARPES by nearest-band RMSE, with the
number of bands per point capped to the true occupied count (flooding is
rejected). PASS < 0.30 eV, PARTIAL 0.30-0.40, FAIL otherwise. A bare KS
submission scores ~0.4-0.6 eV and FAILs by construction.
