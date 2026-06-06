# KS band computation helper for run_qp.py.
# Usage: julia --project=<jenv> _ks_bands.jl <config.json> <grid.csv> <out_ks.csv>
# Emits per grid point: point_id, t, nocc, and the occupied bands plus the
# single lowest unoccupied band (energies in eV relative to E_F), ascending.
using DFTK, PseudoPotentialData, JSON, LinearAlgebra, DelimitedFiles

cfgpath = ARGS[1]; gridpath = ARGS[2]; outpath = ARGS[3]
cfgdir  = dirname(abspath(cfgpath))
cfg = JSON.parsefile(cfgpath)
el  = cfg["element"]; st = cfg["structure"]

# lattice: the three primitive vectors are COLUMNS of `lattice`
latvecs = st["lattice_vectors_bohr"]
lattice = zeros(3,3)
for i in 1:3, j in 1:3
    lattice[j,i] = Float64(latvecs[i][j])   # vector i -> column i
end
positions = [Float64.(p) for p in st["atom_positions_frac"]]   # every atom

dft   = cfg["dft"]
psp   = load_psp(PseudoFamily(dft["pseudopotential_family"]), Symbol(el))
atoms = [ElementPsp(Symbol(el), psp) for _ in positions]
Ecut  = Float64(dft["ecut_Ha"]); kgrid = Int.(dft["kgrid"]); smear = Float64(dft["smearing_Ha"])

model = model_DFT(lattice, atoms, positions;
    functionals=[:lda_x, :lda_c_pw],
    temperature=smear, smearing=DFTK.Smearing.FermiDirac())
basis  = PlaneWaveBasis(model; Ecut, kgrid)
scfres = self_consistent_field(basis; tol=1e-8)
εF     = scfres.εF
zval   = Int(dft["z_valence"]); natoms = length(positions)

# grid points: k_frac = t * endpoint_frac
endpt = Float64.(st["path"]["endpoint_frac"])
grid  = readdlm(gridpath, ',', skipstart=1)
pids  = Int.(grid[:,1]); ts = Float64.(grid[:,2])

# enough bands so occupied + first unoccupied is always available
nb = zval*natoms + 6
b  = compute_bands(scfres, DFTK.ExplicitKpoints([t .* endpt for t in ts]); n_bands=nb)

Ha2eV = 27.211386245988
open(outpath, "w") do io
    println(io, "point_id,t,bands_eV")
    for (i, t) in enumerate(ts)
        ev   = (b.eigenvalues[i] .- εF) .* Ha2eV           # ascending, rel. E_F
        occ  = findall(e -> e < 0.0, ev)
        nocc = isempty(occ) ? 0 : maximum(occ)
        ntake = min(nocc + 1, length(ev))                  # occupied + first unocc
        bands = ev[1:ntake]
        println(io, pids[i], ",", t, ",", join(bands, ";"))
    end
end
