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
