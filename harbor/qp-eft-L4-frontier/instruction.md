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
(`element_config.json` and `grid.csv` only -- no atomic data)
are placed in one directory; your code is given the config and grid paths and
should read the rest relative to the config.

---

# Level 4 — open frontier

Nothing is given but the problem: no formula, no structural ansatz, no prescribed approximations, and no atomic data — you compute your own atomic inputs and choose and justify every truncation. Read `SETUP.md`; its Stage 1 (close the theory on Li in closed form against a theoretical anchor, with a complete contraction enumeration) is mandatory and comes first. The published leading-order treatment is the baseline to match or beat; surpassing it on the concealed metals means physics beyond the published treatment. Document everything in `method.md` (consistency ledger required).

Public development elements: `Li/`, `Na/`, `Al/` — Na and Al with `element_config.json`, `grid.csv`, `arpes_reference.csv`; Li with config + grid only (its self-check anchor is theoretical — see SETUP.md Stage 1, which is mandatory and comes FIRST). NO atomic data files anywhere: whatever atomic inputs your derivation needs, you compute yourself.

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

# Physical setup (Level 4 — open frontier)

## The puzzle and the goal

For the simple metals, Kohn-Sham LDA overestimates the occupied bandwidth seen in
ARPES by 20-35% (alkali) down to a few % (Al). Ordinary DFT functionals (LDA, GGA,
hybrids, G0W0) do not fix it; only a full many-body treatment (eDMFT) reproduces the
measured bandwidths. The missing physics is a DYNAMICAL effect of the core
electrons.

This level is the OPEN rung: nothing about the answer's structure is given — no
propagator ansatz, no decomposition roadmap, no prescribed approximations, and NO
atomic data. You choose the reference, the channels, and every truncation; each
must be declared, controlled, and justified. A published leading-order frozen-core
treatment of this problem exists and is the calibration baseline: through this
exact pinned pipeline it reaches nearest-band RMSE 0.078 eV (Na) and 0.248 eV (Al)
on the public development elements. Matching that level on the concealed metals
means you reconstructed the leading-order physics; beating it there means physics
beyond the published treatment. Both are in scope — the second is the point of
this level.

## Starting point

The rigorous starting point is the ALL-ELECTRON action — valence AND core electrons
in the lattice potential, with the full Coulomb interaction,

    L = int_{r,sigma} bar_psi_{r,sigma} [ d_tau - nabla^2/2 + V_Lat - mu ] psi_{r,sigma}
        + (1/2) int_{r,r'} ( bar_psi_{r,sigma} bar_psi_{r',sigma'}
                              psi_{r',sigma'} psi_{r,sigma} ) / |r - r'| .

Integrate the core fields out of the path integral to obtain an effective action
for the valence electrons alone, and DERIVE — do not posit — the correction it
implies for the quasiparticle band energies. One bookkeeping fact you must
respect: the pinned KS step already contains a static frozen-core pseudopotential
(`V_PSP`), so the static part of whatever your integrate-out yields is already
counted — your correction must consist only of what static pseudopotentials omit,
and re-adding any static piece double-counts the core.

## Stage 1 (mandatory): close the theory on Li before any heavy element

Lithium is provided as a third development element (`Li/`: config + grid; no
ARPES exists for Li). Its core is a single hydrogenic 1s² shell, so EVERY
quantity in your derivation has a closed form. Before generalizing to Na/Al:

1. Derive the full correction on Li ANALYTICALLY: explicit expressions for the
   core orbital, its orbital energy, every Coulomb integral, every core
   excitation energy, the coupling vertex, and z_core at the band bottom. The
   two-electron core makes the enumeration of ALL second-order contractions of
   the core-valence interaction finite — enumerate them COMPLETELY, state what
   each contributes, and carry the complete set forward.
2. Validate against the theoretical anchor: for Li the literature gives a
   Γ-point occupied depth of LDA 3.48 eV vs eDMFT 2.60 eV (no ARPES exists) —
   the many-body narrowing implies z_Γ ≈ 0.75, the LARGEST of the simple
   metals. Your pinned KS setup reproduces the LDA depth to ~0.1 eV; your
   derived correction must reproduce the ~25% narrowing scale from first
   principles. If your derived channels give far less on Li, a contraction is
   missing from your enumeration — return to step 1; do not proceed to heavier
   elements on the strength of Na/Al agreement alone (Li is deliberately the
   element where incomplete channel sets fail loudest).
3. Only then generalize: every step that is closed-form on Li but approximate
   for multi-shell cores (orbital choice, excitation energies, screening) is
   an entry in the approximation ledger with an error estimate.

## What you have — and what you must build yourself

`element_config.json` (lattice, atom positions, `Z_nuclear`, `dft.z_valence`, the
pinned DFT settings), `grid.csv`, and the public ARPES for Na and Al; for Li,
config + grid and the theoretical anchor above (no ARPES). There is NO
atomic core data in this packet: whatever atomic inputs your derivation needs —
core orbitals, core potentials, excitation energies — you compute yourself (e.g.
with your own radial atomic solver; the core occupation follows from `Z_nuclear`
and `z_valence`). Choose the level of atomic theory and justify it. Any core
excitation energy you use must be defined by an explicit expression in your own
reference — state it, do not treat it as a black box.

## Required consistency ledger (in method.md)

1. An approximation ledger: every truncation (reference choice, channel
   selection, closure/spectral treatment, frozen vs relaxed core orbitals,
   screening) declared, with a controlled error estimate.
2. The explicit expression for every core excitation energy you use.
3. For each core channel: the diagonal matrix element of your derived coupling
   against that channel's orbital, and why its value is consistent with your
   energy bookkeeping.
4. If you claim accuracy beyond the published leading order, attribute the gain
   to the specific physical term that produced it.
5. The Li closed-form list: analytic expressions for the Li core orbital
   energy, self-Coulomb integral, excitation energy, coupling vertex, and
   z_Γ — with the numerical comparison against the LDA/eDMFT anchor.

## Rules

Parameter-free: no fitted constant. If you need to tune a constant to match the
public elements, your derivation is incomplete — re-derive, do not fit. The same
code path runs on every element; it is graded on held-out metals whose core
structure differs from the public Na/Al, so matching Na/Al is necessary but not
sufficient.

Your derivation must remain valid across core-shell types — shells both tighter
and more diffuse than the public elements' (e.g. n = 3 s-channels of much larger
radial extent): the shell- and Z-dependence of the coupling must come out of the
derivation, not be calibrated on the tight cores of the development set.

Derive first, evaluate second. Do not use public-element magnitudes to SELECT
between candidate couplings: entire families of incorrect vertices are nearly
degenerate on Na/Al, so numerical agreement there can never substitute for the
derivation. Conversely, a correctly derived coupling can appear wrong by a large
factor if your Bloch-state contraction conventions are inconsistent — if the
magnitude looks off, audit the conventions end-to-end against the public ARPES
before abandoning the derivation.
