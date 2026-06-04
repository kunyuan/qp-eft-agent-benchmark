# Task 01: DFTK Kohn-Sham Baseline

Goal: compute the public Na and Al Kohn-Sham occupied bands using DFTK.

Use:

- LDA
- GTH pseudopotentials
- `ecut_Ha` from `element_config.json`
- Fermi-Dirac smearing `smearing_Ha` from `element_config.json`
- public grid files in `data/public/<element>/grid.csv`

Produce:

```text
ks_bands.csv
```

with columns:

```text
element,point_id,t,band_index,E_KS_minus_EF_eV
```

Write a short note with:

- DFTK version if available;
- pseudopotential family;
- cutoff and smearing;
- the deepest occupied energy at Gamma or the nearest available Gamma point.

This task is a baseline. Do not apply the quasiparticle correction yet.

