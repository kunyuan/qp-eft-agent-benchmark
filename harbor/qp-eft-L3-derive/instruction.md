# Working environment

You are in a container with Julia + DFTK 0.7.25 + PseudoPotentialData and Python
(numpy/scipy) preinstalled; the pinned Julia project is at `$JULIA_PROJECT`.
Public development data is in `./packet/Na/` and `./packet/Al/`. Write your
solution as **`run_qp.py` in the working directory (`/app`)**. The verifier will
invoke it as

```
python run_qp.py --element-config <config.json> --grid <grid.csv> --out <out.csv>
```

on concealed held-out metals (not Na/Al). Each hidden element's data files
(`element_config.json`, `grid.csv`, `core_model.json`, and this level's data)
are placed in one directory; your code is given the config and grid paths and
should read the rest relative to the config.

---

# Level 3 — derive the correction

No formula is given. From the physical setup and the atomic core data, derive the leading frozen-core quasiparticle correction and implement it. Read `SETUP.md`. This is the frontier rung: it amounts to reconstructing the paper's derivation from the stated mechanism and the provided ingredients — expect it to be hard, and document your reasoning in `method.md`.

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


---

# Physical setup (Level 3)

For these simple metals the Kohn-Sham occupied bandwidth overestimates ARPES by
20-35% (alkali) down to a few % (Al). DFT gives no reason KS eigenvalues should
be quasiparticle energies — yet they nearly are. Your task is to find, derive,
and implement the leading parameter-free correction, then predict band energies.

What is missing is dynamical. A frozen-core (large-core) pseudopotential removes
the core electrons; its STATIC core-valence physics is already included. The
omitted piece is the dynamical response of virtual core excitations. Two facts
make the leading correction controlled:

1. Scale separation: core excitation energies `DeltaE_c` (provided per core
   s-channel in `core_model.json`) are several Hartree, far above the valence
   Fermi energy — they can be integrated out as high-energy modes.
2. The interacting uniform electron gas at metallic density is approximately
   Galilean invariant over the occupied Fermi ball, so the static tree-level
   Hamiltonian coincides with the Kohn-Sham Hamiltonian; only the quasiparticle
   pole is rescaled by the frozen-core dynamics.

You are given, per core s-channel c: the all-electron radial core orbital
`u_c(r)` and its single-orbital Hartree potential `V_H_c(r)`
(`atomic_core_<c>.csv`), and the core excitation energy `DeltaE_c`. From your
DFTK calculation you have the KS eigenvalues and the plane-wave coefficients of
each Bloch state. Derive the leading post-SCF quasiparticle correction and
implement it. It must be parameter-free and use the same code path for every
element.
