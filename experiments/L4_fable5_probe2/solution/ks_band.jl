# Pinned Kohn-Sham step: LDA / GTH largecore / settings from element_config.json.
# Usage: julia --project=<env> ks_band.jl <element_config.json> <grid.csv> <out.json>
#
# Outputs JSON with, per grid k-point: KS eigenvalues (Ha), Fermi level (Ha),
# integer G-vectors and complex plane-wave coefficients of each band, plus the
# reciprocal lattice (columns b1,b2,b3, in 1/bohr).

using DFTK
using PseudoPotentialData
using JSON
using LinearAlgebra

function main()
    config_path, grid_path, out_path = ARGS[1], ARGS[2], ARGS[3]
    cfg = JSON.parsefile(config_path)

    element = String(cfg["element"])
    lat = cfg["structure"]["lattice_vectors_bohr"]
    # a1,a2,a3 as COLUMNS of the lattice matrix
    lattice = hcat([Float64.(v) for v in lat]...)
    positions = [Float64.(p) for p in cfg["structure"]["atom_positions_frac"]]
    natoms = length(positions)

    dft = cfg["dft"]
    family = PseudoFamily(String(dft["pseudopotential_family"]))
    psp = load_psp(family, Symbol(element))
    atoms = fill(ElementPsp(Symbol(element), psp), natoms)

    Ecut = Float64(dft["ecut_Ha"])
    kgrid = Int.(dft["kgrid"])
    temperature = Float64(dft["smearing_Ha"])
    zval = Int(dft["z_valence"])

    model = model_DFT(lattice, atoms, positions;
                      functionals=LDA(),
                      temperature, smearing=Smearing.FermiDirac())
    basis = PlaneWaveBasis(model; Ecut, kgrid=MonkhorstPack(kgrid))
    scfres = self_consistent_field(basis; tol=1e-10)
    εF = scfres.εF

    # explicit k-points from the grid
    grid_lines = readlines(grid_path)
    hdr = split(strip(grid_lines[1]), ',')
    @assert hdr[1] == "point_id" && hdr[2] == "t"
    point_ids = Int[]
    ts = Float64[]
    for ln in grid_lines[2:end]
        isempty(strip(ln)) && continue
        f = split(strip(ln), ',')
        push!(point_ids, parse(Int, f[1]))
        push!(ts, parse(Float64, f[2]))
    end
    endpoint = Float64.(cfg["structure"]["path"]["endpoint_frac"])
    kcoords = [t * endpoint for t in ts]

    nb = ceil(Int, zval * natoms / 2) + 6
    band_data = compute_bands(scfres, ExplicitKpoints(kcoords);
                              n_bands=nb, tol=1e-10)

    bs_basis = band_data.basis
    out = Dict{String,Any}()
    out["element"] = element
    out["EF_Ha"] = εF
    out["recip_lattice"] = [collect(c) for c in eachcol(model.recip_lattice)]
    out["lattice"] = [collect(c) for c in eachcol(model.lattice)]
    out["unit_cell_volume"] = model.unit_cell_volume
    out["atom_positions_frac"] = positions
    out["z_valence"] = zval

    kpts_out = Any[]
    for (ik, kpt) in enumerate(bs_basis.kpoints)
        Gs = G_vectors(bs_basis, kpt)
        ψk = band_data.ψ[ik]          # (n_G, n_bands)
        evs = band_data.eigenvalues[ik]
        entry = Dict{String,Any}()
        entry["point_id"] = point_ids[ik]
        entry["t"] = ts[ik]
        entry["kfrac"] = collect(kpt.coordinate)
        entry["G"] = [collect(g) for g in Gs]
        bands = Any[]
        for n = 1:size(ψk, 2)
            push!(bands, Dict("eps_Ha" => evs[n],
                              "c_re" => real.(ψk[:, n]),
                              "c_im" => imag.(ψk[:, n])))
        end
        entry["bands"] = bands
        push!(kpts_out, entry)
    end
    out["kpoints"] = kpts_out

    open(out_path, "w") do io
        JSON.print(io, out)
    end
    println("WROTE ", out_path, "  EF_Ha=", εF)
end

main()
