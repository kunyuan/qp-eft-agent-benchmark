# DFTK Notes

The baseline calculation should use:

- DFTK
- LDA exchange-correlation
- GTH pseudopotentials
- plane-wave cutoff from `element_config.json` (`ecut_Ha`)
- Fermi-Dirac smearing from `element_config.json` (`smearing_Ha`)

The public configs use:

| Element | Structure | Path | Valence | Ecut (Ha) |
| --- | --- | --- | --- | --- |
| Na | bcc | Gamma-N | 1 | 30 |
| Al | fcc | Gamma-X | 3 | 30 |

Use the exact lattice constants in `element_config.json`. The grid file gives
the requested fractional coordinate `t` along the path. Your code should map
`t=0` to Gamma and `t=1` to the named zone-boundary point.

For the submitted runner, it is acceptable to call Julia/DFTK from Python, to
write the whole runner in Julia and wrap it with Python, or to use precomputed
Kohn-Sham outputs generated during your development. The hidden evaluator only
requires the command-line interface documented in `README.md`.

If the benchmark runner is offline, DFTK must already be available in the Julia
environment. The included `environment/Project.toml` documents the intended
Julia dependency.

