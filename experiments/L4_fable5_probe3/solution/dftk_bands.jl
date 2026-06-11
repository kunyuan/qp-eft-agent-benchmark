# Pinned-DFTK runner for the quasiparticle pipeline (element-blind).
# Runs the pinned LDA/GTH SCF, computes bands at the explicit k-points
# k_frac = t * endpoint_frac, and dumps eigenvalues, Fermi level, cell data and
# plane-wave coefficients (cartesian k+G) for the core-coupling evaluation.
#
# Usage: julia --project=<pinned env> dftk_bands.jl config.json tlist.txt outdir nbands
using DFTK, PseudoPotentialData, JSON, LinearAlgebra, Printf

cfgpath, tpath, outdir = ARGS[1], ARGS[2], ARGS[3]
nbands = parse(Int, ARGS[4])
cfg = JSON.parsefile(cfgpath)
mkpath(outdir)

elem = Symbol(cfg["element"])
lattice = hcat([Float64.(v) for v in cfg["structure"]["lattice_vectors_bohr"]]...)
positions = [Float64.(p) for p in cfg["structure"]["atom_positions_frac"]]
endpoint = Float64.(cfg["structure"]["path"]["endpoint_frac"])
dft = cfg["dft"]
psp = load_psp(PseudoFamily(dft["pseudopotential_family"]), elem)
atoms = [ElementPsp(elem, psp) for _ in positions]
model = model_DFT(lattice, atoms, positions; functionals=LDA(),
                  temperature=Float64(dft["smearing_Ha"]),
                  smearing=Smearing.FermiDirac())
basis = PlaneWaveBasis(model; Ecut=Float64(dft["ecut_Ha"]), kgrid=Int.(dft["kgrid"]))
scfres = self_consistent_field(basis; tol=1e-8)

ts = [parse(Float64, strip(l)) for l in readlines(tpath) if !isempty(strip(l))]
kcoords = [t * endpoint for t in ts]
bands = compute_bands(scfres, ExplicitKpoints(kcoords); n_bands=nbands, tol=1e-9)

open(joinpath(outdir, "ef.txt"), "w") do io
    println(io, scfres.εF)
end
recip = model.recip_lattice
open(joinpath(outdir, "bands.csv"), "w") do io
    println(io, "ik,t,n,eps_Ha")
    for (ik, t) in enumerate(ts), n in 1:nbands
        @printf(io, "%d,%.8f,%d,%.10f\n", ik, t, n, bands.eigenvalues[ik][n])
    end
end
open(joinpath(outdir, "psi.csv"), "w") do io
    println(io, "ik,n,kg,w")   # |k+G| (bohr^-1) and |c_G|^2
    for (ik, kfrac) in enumerate(kcoords)
        kpt = bands.basis.kpoints[ik]
        Gs = G_vectors(bands.basis, kpt)
        kgmag = [norm(recip * (kfrac + G)) for G in Gs]
        psik = bands.ψ[ik]
        for n in 1:nbands
            c = psik[:, n]
            for (ig, kg) in enumerate(kgmag)
                w = abs2(c[ig])
                w > 1e-12 && @printf(io, "%d,%d,%.8f,%.10e\n", ik, n, kg, w)
            end
        end
    end
end
open(joinpath(outdir, "cell.txt"), "w") do io
    println(io, abs(det(lattice)))
    println(io, length(positions))
end
println("DFTK_BANDS_OK")
