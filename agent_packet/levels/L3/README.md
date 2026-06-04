# Level 3 — derive the correction

No formula is given. From the physical setup and the atomic core data, derive the leading frozen-core quasiparticle correction and implement it. Read `SETUP.md`.

Public development elements: `Na/`, `Al/` (each with `element_config.json`, `grid.csv`, `arpes_reference.csv`, and the level's data files).

## Submission

Submit a directory containing `run_qp.py` with this exact interface:

```bash
python run_qp.py --element-config element_config.json --grid grid.csv --out qp_bands.csv
```

Output CSV columns: `element,point_id,t,E_pred_eV` — one row per occupied
quasiparticle band at each grid point, energy relative to the Fermi level.
`t` in `grid.csv` is the fractional path coordinate: `k_frac = t * endpoint_frac`
(both in `element_config.json`).

## Pinned DFT setup (use exactly; do not "improve" it)

LDA / GTH `cp2k.nc.sr.lda.v0_1.largecore.gth` / `Ecut_Ha`, `kgrid`,
`smearing_Ha` from the config / Fermi-Dirac smearing. The Kohn-Sham band is
the starting point — a converged KS band alone is NOT a quasiparticle
prediction for this benchmark.

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
