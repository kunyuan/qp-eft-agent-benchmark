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

Build the crystal from the config: `structure.lattice_vectors_bohr` are the
three primitive lattice vectors a1,a2,a3 in Bohr (in DFTK use them as the COLUMNS
of the `lattice` matrix), and `structure.atom_positions_frac` are the fractional
atom coordinates — there may be MORE THAN ONE atom (e.g. hcp has two), so do not
assume a single-atom cell.

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

At each grid point, with the bands in ascending energy order, emit **every
occupied band (`E_KS < E_F`) plus the single lowest unoccupied band**, then stop.
(No energy-margin threshold: including the first unoccupied band captures a Fermi
crossing — where ARPES can resolve a band just below E_F that DFT places just
above — and is robust to the exact position of E_F at the edge.) Report
`E_pred_eV = E_QP - E_F` (eV); it may be slightly positive for that first
unoccupied band. Compute enough bands (a few above `z_valence * n_atoms`) so the
occupied set plus one is available. The evaluator caps predictions to the true
band count per point, so do not emit extra bands.

## Rules

- Parameter-free: no fitting to ARPES, no per-element tuning, no hardcoded
  output values, no per-element `if` branches. The SAME code path runs on the
  public elements and on concealed held-out metals.
- Develop only against the public elements here (Na, Al). Hidden metals are
  scored by the evaluator and are not revealed.
- **No network access in `run_qp.py`.** Compute everything from the provided
  inputs + DFTK, which runs fully offline (the pinned Julia environment and the
  GTH pseudopotentials are already installed locally). Do not download data, clone
  repositories, query web services, or otherwise reach the network at run time — a
  submission whose prediction depends on a network call is invalid.
- **No external lookup at all — solve from this packet only.** Treat yourself as
  having no internet: do **not** use web search for anything — not the underlying
  research paper / its derivation / the answers, and not library or API
  documentation. Everything you need is here plus the tools installed locally.
  DFTK is installed in the pinned project; discover its exact API by introspecting
  the package in Julia (`names(DFTK; all=true)`, `?<name>`, `methods(<fn>)`, or
  reading the installed source under the depot) — not from the web.

## Accuracy target & scoring

Your predictions are compared to held-out ARPES on the concealed metals by
nearest-band RMSE. You must cover **every** grid point with **exactly** the bands
specified above (occupied + first unoccupied) — flooded, sparse, or wrong-count
submissions are rejected. At evaluation your runner is given inputs WITHOUT the
held-out ARPES, so the prediction must come from the physics.

How accurate is good enough? Bare Kohn–Sham overestimates the occupied bandwidth
by 20–35%; a correct parameter-free frozen-core correction removes essentially all
of that, reaching agreement with experiment at the level of the state-of-the-art
many-body method eDMFT. For scale, Γ-point occupied-band depths (eV below E_F):

| element | LDA (Mandal 2022) | eDMFT (Mandal 2022) | experiment |
|---------|-------------------|---------------------|------------|
| Li | 3.48 | 2.60 | — (no ARPES) |
| Na | 3.30 | 2.84 | 2.65–2.78 |
| Ca | 3.98 | 3.24 | 3.30 |

The LDA column is a literature reference; your pinned DFTK setup reproduces it to
within ~0.1 eV — **use your own computed KS, do not tune your setup to match these
numbers.** Aim for agreement with experiment/eDMFT at the ~0.1 eV level (eDMFT
itself differs from experiment by ~0.1 eV, so this is the realistic target, not an
exact match). Self-check your method against these depths and against the full
ARPES bands provided for Na and Al before it is run on the concealed metals.


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
