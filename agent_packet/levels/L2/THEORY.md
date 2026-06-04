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
