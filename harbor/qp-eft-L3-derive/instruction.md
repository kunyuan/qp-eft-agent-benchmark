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

## The puzzle and the goal

For the simple metals, Kohn-Sham LDA overestimates the occupied bandwidth seen in
ARPES by 20-35% (alkali) down to a few % (Al). Ordinary DFT functionals (LDA, GGA,
hybrids, G0W0) do not fix it; only a full many-body treatment (eDMFT) reproduces the
measured bandwidths. The missing physics is therefore a DYNAMICAL effect of the core
electrons. Your goal: DERIVE, from first principles, the parameter-free correction
that captures it — reaching eDMFT / experiment-level accuracy — and predict the
band energies.

## The quasiparticle propagator (the result you are reproducing)

This effective theory rests on (i) the scale separation between the core excitation
energies and the valence Fermi energy (which lets the core be integrated out), and
(ii) the locality of the interacting uniform electron gas at metallic densities
(established by high-order diagrammatic Monte Carlo), which dictates how many-body
effects are absorbed into renormalized parameters. The tree-level quasiparticle
propagator of the EFT is

    G_{nu,k}^{-1}(i omega)  ~=  (z_val)^{-1} [ (z_core_{nu,k})^{-1} i omega
                                              - ( eps_KS_{nu,k} - eps_F ) ]      (1)

where the dispersion is governed by the eigenvalues `eps_KS_{nu,k}` of the
Kohn-Sham Hamiltonian

    H_KS  =  -1/2 nabla^2  +  V_PSP  +  V_H[n]  +  V_xc[n]      (V_PSP: GTH large-core),

`z_val` is the valence quasiparticle residue (an overall prefactor that does NOT
shift the pole), and `z_core_{nu,k}` is the frozen-core renormalization, which
rescales the frequency relative to the dispersion — an effective time dilation from
virtual core excitations. Locating the pole gives the energy ARPES measures,

    eps_QP_{nu,k}  ~=  eps_F  +  z_core_{nu,k} ( eps_KS_{nu,k} - eps_F ).      (2)

Your DFTK run gives `eps_KS_{nu,k}` and `eps_F`. What is NOT given is
`z_core_{nu,k}` — deriving its closed form is the task.

## How to derive the effective valence action

The rigorous starting point is the ALL-ELECTRON action — valence AND core electrons
in the lattice potential, with the full Coulomb interaction,

    L = int_{r,sigma} bar_psi_{r,sigma} [ d_tau - nabla^2/2 + V_Lat - mu ] psi_{r,sigma}
        + (1/2) int_{r,r'} ( bar_psi_{r,sigma} bar_psi_{r',sigma'}
                              psi_{r',sigma'} psi_{r,sigma} ) / |r - r'| .

INTEGRATE THE CORE FIELDS OUT of the path integral to obtain an effective action for
the valence electrons alone. (A clean route: pass to a reference in which the valence
orbitals are frozen — pushed to high energy by a hybridization — so the core is
solved in that reference and the valence then propagates in the potential the core
leaves behind.) The result is a valence propagator

    g_v^{-1}  =  g_0^{-1}  +  delta_V_pp ,        g_0^{-1} = i omega - H_0,
                                                  H_0 = -nabla^2/2 + V_Lat - mu,

where `delta_V_pp` is the pseudopotential induced by the core. It has
- a STATIC (energy-independent) part: the conventional frozen-core pseudopotential
  `V_PSP` already inside `H_KS` — this sets the dispersion `eps_KS` in Eq. (1);
- a DYNAMICAL (frequency-dependent) part: this is the `z_core` renormalization in
  Eq. (1), and it is the leading thing static pseudopotentials omit.

Re-adding the static part would double-count the core; the correction is the
dynamical part only. Derive it — and hence `z_core_{nu,k}` and, through Eq. (2), the
band energies (relative to E_F) and the bandwidth.

Do this as a real derivation, not a model: start from the many-electron action and
integrate the core fields out STEP BY STEP to obtain `δV_pp` — and in particular its
dynamical (frequency-dependent) term — from the valence–core Coulomb interaction
itself. The coupling that enters `z_core` is whatever this integral yields; do not
posit a model self-energy or fix the coupling by dimensional analysis.

Be careful with EXCHANGE: valence and core electrons are antisymmetric, so `δV_pp`
carries a Fock/exchange contribution alongside the direct (Hartree) one, and the
antisymmetry shapes the dynamical coupling — a purely direct density-response
picture misses it.

One simplification is controlled: the core excitation energies `DeltaE_c` are several
Hartree, far above the valence Fermi energy, so treat the core excitations by closure
at the single scale `DeltaE_c`, and — since they dominate — keep only the core
s-shell (monopole) channels.

Do not posit the form of `z_core_{nu,k}`: what it is built from out of the core data,
and how it depends on the Bloch state, are for you to derive.

## What you have

Per core s-channel `c`: the all-electron radial core orbital `u_c(r)`
(`atomic_core_<c>.csv`, normalized int u_c^2 dr = 1) and its excitation energy
`DeltaE_c` (`core_model.json`). `DeltaE_c` is the full many-electron excitation
energy of the interacting core — the total-energy difference
E_core(hole in c) − E_core(ground state) — NOT a bare orbital-eigenvalue
difference; keep your coupling and your energy denominators consistent with this
definition. From your DFTK run: the KS eigenvalues, E_F, and the
plane-wave coefficients `c_nk(G)` of each Bloch state (`|k+G|` in Bohr^-1). Anything
else your derivation needs — potentials, integrals — you build from these.

Parameter-free: no fitted constant (the correct prefactor is exactly 1). If you need
to tune a constant to match the public elements, your derivation is incomplete —
re-derive, do not fit. The same code path runs on every element; it is graded on
held-out metals whose core structure differs from the public Na/Al, so matching
Na/Al is necessary but not sufficient.

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
