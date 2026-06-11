# Pinned KS pipeline dump: SCF + bands at explicit path k-points.
# Usage: julia --project=<env> ks_dump.jl element_config.json grid.csv out.json
#
# Dumps JSON with: EF, cell metadata, SCF-grid eigenvalues (+weights, for the
# corrected-Fermi-level refill) and, per path k-point and band, the eigenvalue
# and the dominant plane-wave coefficients c_nk(G) (integer G + complex c).
using DFTK, PseudoPotentialData, JSON, LinearAlgebra

cfgpath, gridpath, outpath = ARGS[1], ARGS[2], ARGS[3]

cfg = JSON.parsefile(cfgpath)
elem = cfg["element"]
lat = cfg["structure"]["lattice_vectors_bohr"]
lattice = hcat([Float64.(v) for v in lat]...)          # columns a1,a2,a3 (Bohr)
positions = [Float64.(p) for p in cfg["structure"]["atom_positions_frac"]]
endpoint = Float64.(cfg["structure"]["path"]["endpoint_frac"])
dft = cfg["dft"]

psp = load_psp(PseudoFamily(dft["pseudopotential_family"]), Symbol(elem))
atoms = [ElementPsp(Symbol(elem), psp) for _ in positions]
model = model_DFT(lattice, atoms, positions;
                  functionals=LDA(),
                  temperature=dft["smearing_Ha"],
                  smearing=DFTK.Smearing.FermiDirac())
basis = PlaneWaveBasis(model; Ecut=dft["ecut_Ha"], kgrid=Int.(dft["kgrid"]))
scfres = self_consistent_field(basis; tol=1e-10)

# ---- grid ----
ids = Int[]; ts = Float64[]
for line in readlines(gridpath)[2:end]
    isempty(strip(line)) && continue
    parts = split(line, ",")
    push!(ids, parse(Int, parts[1])); push!(ts, parse(Float64, parts[2]))
end
kcoords = [t * endpoint for t in ts]

nval = sum(DFTK.n_elec_valence(a) for a in atoms)
nbands = Int(ceil(nval)) + 6
bands = compute_bands(scfres, DFTK.ExplicitKpoints(kcoords); n_bands=nbands, tol=1e-9)

NKEEP = 48  # dominant PW coefficients kept per state

path_pts = []
for (i, kc) in enumerate(kcoords)
    kpt = bands.basis.kpoints[i]
    Gs = G_vectors(bands.basis, kpt)            # Vec3{Int} list
    ψk = bands.ψ[i]                              # nG × nbands
    bandlist = []
    for n in 1:nbands
        c = ψk[:, n]
        w = abs2.(c)
        order = sortperm(w; rev=true)
        keep = order[1:min(NKEEP, length(order))]
        push!(bandlist, Dict(
            "eps_Ha" => bands.eigenvalues[i][n],
            "G" => [Int.(Gs[j]) for j in keep],
            "c_re" => [real(c[j]) for j in keep],
            "c_im" => [imag(c[j]) for j in keep],
            "w_kept" => sum(w[keep])))
    end
    push!(path_pts, Dict("point_id" => ids[i], "t" => ts[i],
                         "k_frac" => kc, "bands" => bandlist))
end

# ---- SCF k-grid eigenvalues for Fermi refill ----
scf_grid = []
for (ik, kpt) in enumerate(basis.kpoints)
    push!(scf_grid, Dict("kweight" => basis.kweights[ik],
                         "eps_Ha" => scfres.eigenvalues[ik]))
end

out = Dict(
    "element" => elem,
    "EF_Ha" => scfres.εF,
    "n_electrons" => model.n_electrons,
    "volume_bohr3" => model.unit_cell_volume,
    "lattice_cols" => [lattice[:, i] for i in 1:3],
    "recip_cols" => [model.recip_lattice[:, i] for i in 1:3],
    "positions_frac" => positions,
    "n_atoms" => length(positions),
    "temperature_Ha" => dft["smearing_Ha"],
    "nbands" => nbands,
    "path" => path_pts,
    "scf_grid" => scf_grid)

open(outpath, "w") do io
    JSON.print(io, out)
end
println("wrote ", outpath)
