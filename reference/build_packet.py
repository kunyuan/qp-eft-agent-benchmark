#!/usr/bin/env python3
"""Assemble the four-level agent packet + hidden eval sets from the generated
packet_data and the existing grids/ARPES.

Levels (same physics, withhold ladder):
  L1 apply       : full z_core formula given + precomputed f_c(K) tables.
  L2 formfactor  : formula given, NO f_c tables; atomic core data -> compute f_c.
  L3 derive      : NO formula, only physical setup; atomic core data.
  L4 frontier    : NO formula, NO atomic data; agent computes its own atomic
                   inputs and may go beyond the published leading order.

Layout produced:
  agent_packet/levels/L{1,2,3,4}/{Na,Al}/   (public dev: config, grid, arpes, data)
  evaluator/hidden/L{1,2,3,4}/{K,Mg}/        (hidden: config, grid, arpes[scoring], data)
"""
from __future__ import annotations
import json, math, shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PDATA = ROOT / "reference" / "packet_data"

# pinned per-element setup (matches reference/gold_runner.jl ELEMENTS & DESIGN.md)
META = {
    "Na": dict(role="public", Z=11, bravais="bcc", a=8.107, ca=None, zval=1,
               ecut=15.0, kgrid=[8, 8, 8], path=("Gamma", "N", [0, 0, 0.5]),
               core_s=["1s", "2s"], valence="3s1",
               src_grid=ROOT/"evaluator/source_data/Na"),
    "Al": dict(role="public", Z=13, bravais="fcc", a=7.653, ca=None, zval=3,
               ecut=20.0, kgrid=[8, 8, 8], path=("Gamma", "X", [0.5, 0, 0.5]),
               core_s=["1s", "2s"], valence="3s2 3p1",
               src_grid=ROOT/"evaluator/source_data/Al"),
    "K":  dict(role="hidden", Z=19, bravais="bcc", a=9.874, ca=None, zval=1,
               ecut=15.0, kgrid=[8, 8, 8], path=("Gamma", "N", [0, 0, 0.5]),
               core_s=["1s", "2s", "3s"], valence="4s1",
               src_grid=ROOT/"evaluator/source_data/K"),
    "Mg": dict(role="hidden", Z=12, bravais="hcp", a=6.066, ca=1.624, zval=2,
               ecut=20.0, kgrid=[8, 8, 6], path=("Gamma", "A", [0, 0, 0.5]),
               core_s=["1s", "2s"], valence="3s2",
               src_grid=ROOT/"evaluator/source_data/Mg"),
}


def lattice_and_positions(m: dict):
    """Explicit primitive cell (matches the gold's construction exactly).
    Returns (lattice_vectors_bohr, atom_positions_frac):
    lattice_vectors_bohr = [a1, a2, a3] = the COLUMNS of the DFTK lattice matrix."""
    a = m["a"]
    if m["bravais"] == "bcc":
        vecs = [[-a/2, a/2, a/2], [a/2, -a/2, a/2], [a/2, a/2, -a/2]]
        pos = [[0.0, 0.0, 0.0]]
    elif m["bravais"] == "fcc":
        vecs = [[0.0, a/2, a/2], [a/2, 0.0, a/2], [a/2, a/2, 0.0]]
        pos = [[0.0, 0.0, 0.0]]
    elif m["bravais"] == "hcp":
        c = a * m["ca"]
        vecs = [[a, 0.0, 0.0], [a/2, a * math.sqrt(3) / 2, 0.0], [0.0, 0.0, c]]
        pos = [[0.0, 0.0, 0.0], [1/3, 2/3, 1/2]]
    else:
        raise ValueError(f"unsupported bravais {m['bravais']}")
    return vecs, pos


def config_for(el: str, m: dict) -> dict:
    vecs, pos = lattice_and_positions(m)
    struct = {"bravais": m["bravais"], "a_bohr": m["a"],
              "lattice_vectors_bohr": vecs, "atom_positions_frac": pos,
              "path": {"from": m["path"][0], "to": m["path"][1],
                       "endpoint_frac": m["path"][2]}}
    if m["ca"] is not None:
        struct["c_over_a"] = m["ca"]
    return {
        "element": el,
        "Z_nuclear": m["Z"],
        "structure": struct,
        "dft": {"xc": "LDA",
                "pseudopotential_family": "cp2k.nc.sr.lda.v0_1.largecore.gth",
                "z_valence": m["zval"], "ecut_Ha": m["ecut"],
                "kgrid": m["kgrid"], "smearing_Ha": 0.001},
        "frozen_core": {"core_s_channels": m["core_s"], "valence": m["valence"],
                        "core_model_file": "core_model.json"},
        "note": ("Build the cell from structure.lattice_vectors_bohr (the three "
                 "primitive lattice vectors a1,a2,a3 in Bohr; in DFTK use them as "
                 "the COLUMNS of the `lattice` matrix) and structure.atom_positions_frac "
                 "(fractional coordinates, one per atom -- NOT necessarily a single "
                 "atom). t in grid.csv is fractional along Gamma->endpoint: "
                 "k_frac = t * endpoint_frac."),
    }


def _write_core_orbital_only(src_csv: Path, dst_csv: Path):
    """L3: emit atomic_core CSV with only r_bohr,u_c (drop the derivable V_H_c)."""
    rows = src_csv.read_text().splitlines()
    out = ["r_bohr,u_c"]
    for row in rows[1:]:
        if not row.strip():
            continue
        c = row.split(",")
        out.append(f"{c[0]},{c[1]}")
    dst_csv.write_text("\n".join(out) + "\n")


def copy_level_data(el: str, level: int, dst: Path):
    if level == 4:
        # L4 (open frontier): NO atomic data at all. The agent computes its own
        # core orbitals / potentials / excitation energies from Z_nuclear and
        # z_valence (already in element_config.json) with its own atomic solver.
        return
    src = PDATA / el
    full = json.loads((src / "core_model.json").read_text())
    # level-trimmed core model: L1 gives f_c tables (+ DeltaE); L2/L3 give only
    # DeltaE_c + channel identity (f0/J_c are form-factor outputs to be computed).
    chans = []
    for c in full["core_s_channels"]:
        entry = {"channel": c["channel"], "DeltaE_c_Ha": c["DeltaE_c_Ha"]}
        if level == 1:
            entry["fc_table_file"] = f"fc_table_{c['channel']}.csv"
        else:
            entry["atomic_core_file"] = f"atomic_core_{c['channel']}.csv"
        chans.append(entry)
    cm = {"element": el, "core_s_channels": chans,
          "DeltaE_c_note": ("DeltaE_c_Ha is the full interacting-core excitation "
                            "energy E_core(hole in c) - E_core(ground state), a "
                            "many-electron total-energy difference -- NOT a bare "
                            "orbital-eigenvalue difference.")}
    if level == 2:
        # L2 gives the f_c formula explicitly, and that formula uses V_H_c.
        cm["atomic_data_columns"] = {"atomic_core_*.csv": ["r_bohr", "u_c", "V_H_c"]}
        cm["atomic_data_note"] = ("u_c(r)=r*R_c(r) normalized so integral u_c^2 dr = 1; "
                                  "V_H_c is the single-core-orbital Hartree potential.")
    elif level == 3:
        # L3 derives from scratch: hand over only the raw core orbital u_c. The
        # agent must build any potential it needs (e.g. V_H_c) from u_c itself --
        # which V_H_c is the coupling is part of what the derivation has to find.
        cm["atomic_data_columns"] = {"atomic_core_*.csv": ["r_bohr", "u_c"]}
        cm["atomic_data_note"] = "u_c(r)=r*R_c(r) normalized so integral u_c^2 dr = 1."
    (dst / "core_model.json").write_text(json.dumps(cm, indent=2) + "\n")
    for ch in META[el]["core_s"]:
        if level == 1:
            shutil.copy(src / f"fc_table_{ch}.csv", dst / f"fc_table_{ch}.csv")
        elif level == 3:
            _write_core_orbital_only(src / f"atomic_core_{ch}.csv",
                                     dst / f"atomic_core_{ch}.csv")
        else:
            shutil.copy(src / f"atomic_core_{ch}.csv", dst / f"atomic_core_{ch}.csv")


ENVIRONMENT = """\
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
"""
SUBMISSION = ENVIRONMENT

THEORY_L1 = """\
# Theory (Level 1)

The leading frozen-core quasiparticle correction is post-SCF: compute the KS
band, then compress occupied energies toward the Fermi level by a state-
dependent factor.

```
E_QP(n,k) - E_F = z_core(n,k) * (E_KS(n,k) - E_F)
z_core(n,k)     = 1 / (1 + sum_c |F_c(n,k)|^2 / DeltaE_c^2)
F_c(n,k)        = sum_G c_nk(G) * f_c(|k+G|)
```

- `c_nk(G)` are the plane-wave coefficients of the KS Bloch state |n,k>
  (read them from your DFTK calculation).
- `f_c(K)` is the core form factor for channel c — **given** in
  `fc_table_<c>.csv` (columns `K_bohr_inv,f_c`); interpolate.
- `DeltaE_c` (Ha) is in `core_model.json`.
- `|k+G|` is the Cartesian length in Bohr^-1 (use the reciprocal lattice).

Do not add the static core self-energy again — it is already in the
pseudopotential. Only the frequency-dependent piece (the z_core factor) is new.
"""

THEORY_L2 = """\
# Theory (Level 2)

Same correction as Level 1:

```
E_QP(n,k) - E_F = z_core(n,k) * (E_KS(n,k) - E_F)
z_core(n,k)     = 1 / (1 + sum_c |F_c(n,k)|^2 / DeltaE_c^2)
F_c(n,k)        = sum_G c_nk(G) * f_c(|k+G|)
```

but the core form factor `f_c(K)` is **not** given. Compute it from the atomic
core data (`atomic_core_<c>.csv`: `r_bohr,u_c,V_H_c`):

```
J_c     = integral u_c(r)^2 V_H_c(r) dr
f_c(K)  = sqrt(4*pi)/K * integral u_c(r) [V_H_c(r) - J_c] sin(K r) dr     (K>0)
f_c(0)  = sqrt(4*pi)   * integral u_c(r) [V_H_c(r) - J_c] r dr
```

`DeltaE_c` (Ha) is in `core_model.json`. `c_nk(G)` come from your DFTK run;
`|k+G|` is the Cartesian length in Bohr^-1.
"""

SETUP_L3 = """\
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
"""

SETUP_L4 = """\
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

## What you have — and what you must build yourself

`element_config.json` (lattice, atom positions, `Z_nuclear`, `dft.z_valence`, the
pinned DFT settings), `grid.csv`, and the public ARPES for Na and Al. There is NO
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
"""

LEVEL_DOCS = {
    1: ("Level 1 — apply the given correction",
        "The full formula and the core form factors are given. Wire up DFTK with "
        "the pinned settings, extract the Bloch coefficients, assemble `z_core`, "
        "and predict. Read `THEORY.md`.", "THEORY.md", THEORY_L1),
    2: ("Level 2 — compute the form factor",
        "The formula is given, but you must compute the core form factor `f_c(K)` "
        "yourself from the atomic core data. Read `THEORY.md`.", "THEORY.md", THEORY_L2),
    3: ("Level 3 — derive the correction",
        "No formula is given. From the physical setup and the atomic core data, "
        "derive the leading frozen-core quasiparticle correction and implement it. "
        "Read `SETUP.md`. This is the frontier rung: it amounts to reconstructing "
        "the paper's derivation from the stated mechanism and the provided "
        "ingredients — expect it to be hard, and document your reasoning in "
        "`method.md`.", "SETUP.md", SETUP_L3),
    4: ("Level 4 — open frontier",
        "Nothing is given but the problem: no formula, no structural ansatz, no "
        "prescribed approximations, and no atomic data — you compute your own "
        "atomic inputs and choose and justify every truncation. Read `SETUP.md`. "
        "The published leading-order treatment is the baseline to match or beat; "
        "surpassing it on the concealed metals means physics beyond the published "
        "treatment. Document everything in `method.md` (consistency ledger "
        "required).", "SETUP.md", SETUP_L4),
}


def write_docs():
    base = ROOT / "agent_packet" / "levels"
    for level, (title, blurb, docname, doctext) in LEVEL_DOCS.items():
        d = base / f"L{level}"
        (d / docname).write_text(doctext)
        data_phrase = ("and NO atomic data files — whatever atomic inputs your "
                       "derivation needs, you compute yourself" if level == 4
                       else "and the level's data files")
        (d / "README.md").write_text(
            f"# {title}\n\n{blurb}\n\nPublic development elements: `Na/`, `Al/` "
            f"(each with `element_config.json`, `grid.csv`, `arpes_reference.csv`, "
            f"{data_phrase}).\n\n{SUBMISSION}")
    (ROOT / "agent_packet" / "README.md").write_text(
        "# Frozen-Core Quasiparticle Band Benchmark — agent packet\n\n"
        "Predict occupied quasiparticle band energies for simple metals. The same\n"
        "physics problem is offered at four difficulty levels (pick one); all are\n"
        "scored identically against held-out ARPES.\n\n"
        "- `levels/L1/` — apply the given `z_core` formula + given form factors.\n"
        "- `levels/L2/` — formula given; compute the form factor from atomic data.\n"
        "- `levels/L3/` — derive the correction from the physical setup.\n"
        "- `levels/L4/` — open frontier: no formula, no atomic data; compute your\n"
        "  own atomic inputs, declare every approximation, beat the published\n"
        "  leading-order baseline if you can.\n\n"
        "Develop against the public elements (Na, Al). The evaluator runs your\n"
        "`run_qp.py` on concealed held-out metals via the same interface.\n\n"
        + SUBMISSION)


def build():
    # clean previous build
    if (ROOT/"agent_packet/levels").exists():
        shutil.rmtree(ROOT/"agent_packet/levels")
    for level in (1, 2, 3, 4):
        for h in ("K", "Mg"):
            hp = ROOT/f"evaluator/hidden/L{level}"/h
            if hp.exists():
                shutil.rmtree(hp)

    for el, m in META.items():
        cfg = config_for(el, m)
        for level in (1, 2, 3, 4):
            if m["role"] == "public":
                dst = ROOT/f"agent_packet/levels/L{level}"/el
            else:
                dst = ROOT/f"evaluator/hidden/L{level}"/el
            dst.mkdir(parents=True, exist_ok=True)
            cfg_lvl = {k: v for k, v in cfg.items() if not (level == 4 and k == "frozen_core")}
            (dst/"element_config.json").write_text(json.dumps(cfg_lvl, indent=2) + "\n")
            shutil.copy(m["src_grid"]/"grid.csv", dst/"grid.csv")
            shutil.copy(m["src_grid"]/"arpes_reference.csv", dst/"arpes_reference.csv")
            copy_level_data(el, level, dst)
    write_docs()
    print("packet assembled")


if __name__ == "__main__":
    build()
