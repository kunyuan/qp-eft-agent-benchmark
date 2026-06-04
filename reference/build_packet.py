#!/usr/bin/env python3
"""Assemble the three-level agent packet + hidden eval sets from the generated
packet_data and the existing grids/ARPES.

Levels (same physics, withhold ladder):
  L1 apply       : full z_core formula given + precomputed f_c(K) tables.
  L2 formfactor  : formula given, NO f_c tables; atomic core data -> compute f_c.
  L3 derive      : NO formula, only physical setup; atomic core data.

Layout produced:
  agent_packet/levels/L{1,2,3}/{Na,Al}/   (public dev: config, grid, arpes, data)
  evaluator/hidden/L{1,2,3}/{K,Mg}/        (hidden: config, grid, arpes[scoring], data)
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


def copy_level_data(el: str, level: int, dst: Path):
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
    cm = {"element": el, "core_s_channels": chans}
    if level >= 2:
        cm["atomic_data_columns"] = {"atomic_core_*.csv": ["r_bohr", "u_c", "V_H_c"]}
        cm["atomic_data_note"] = ("u_c(r)=r*R_c(r) normalized so integral u_c^2 dr = 1; "
                                  "V_H_c is the single-core-orbital Hartree potential.")
    (dst / "core_model.json").write_text(json.dumps(cm, indent=2) + "\n")
    for ch in META[el]["core_s"]:
        if level == 1:
            shutil.copy(src / f"fc_table_{ch}.csv", dst / f"fc_table_{ch}.csv")
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

## Scoring

Predictions are compared to held-out ARPES by nearest-band RMSE. You must cover
**every** grid point with **exactly** the bands specified above (occupied + first
unoccupied) — flooded, sparse, or wrong-count submissions are rejected. PASS <
0.30 eV, PARTIAL 0.30-0.40, FAIL otherwise. A bare KS submission scores ~0.4-0.6
eV and FAILs by construction. (At evaluation your runner is given inputs WITHOUT
the held-out ARPES, so the prediction must come from the physics.)
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
}


def write_docs():
    base = ROOT / "agent_packet" / "levels"
    for level, (title, blurb, docname, doctext) in LEVEL_DOCS.items():
        d = base / f"L{level}"
        (d / docname).write_text(doctext)
        (d / "README.md").write_text(
            f"# {title}\n\n{blurb}\n\nPublic development elements: `Na/`, `Al/` "
            f"(each with `element_config.json`, `grid.csv`, `arpes_reference.csv`, "
            f"and the level's data files).\n\n{SUBMISSION}")
    (ROOT / "agent_packet" / "README.md").write_text(
        "# Frozen-Core Quasiparticle Band Benchmark — agent packet\n\n"
        "Predict occupied quasiparticle band energies for simple metals. The same\n"
        "physics problem is offered at three difficulty levels (pick one); all are\n"
        "scored identically against held-out ARPES.\n\n"
        "- `levels/L1/` — apply the given `z_core` formula + given form factors.\n"
        "- `levels/L2/` — formula given; compute the form factor from atomic data.\n"
        "- `levels/L3/` — derive the correction from the physical setup.\n\n"
        "Develop against the public elements (Na, Al). The evaluator runs your\n"
        "`run_qp.py` on concealed held-out metals via the same interface.\n\n"
        + SUBMISSION)


def build():
    # clean previous build
    if (ROOT/"agent_packet/levels").exists():
        shutil.rmtree(ROOT/"agent_packet/levels")
    for level in (1, 2, 3):
        for h in ("K", "Mg"):
            hp = ROOT/f"evaluator/hidden/L{level}"/h
            if hp.exists():
                shutil.rmtree(hp)

    for el, m in META.items():
        cfg = config_for(el, m)
        for level in (1, 2, 3):
            if m["role"] == "public":
                dst = ROOT/f"agent_packet/levels/L{level}"/el
            else:
                dst = ROOT/f"evaluator/hidden/L{level}"/el
            dst.mkdir(parents=True, exist_ok=True)
            (dst/"element_config.json").write_text(json.dumps(cfg, indent=2) + "\n")
            shutil.copy(m["src_grid"]/"grid.csv", dst/"grid.csv")
            shutil.copy(m["src_grid"]/"arpes_reference.csv", dst/"arpes_reference.csv")
            copy_level_data(el, level, dst)
    write_docs()
    print("packet assembled")


if __name__ == "__main__":
    build()
