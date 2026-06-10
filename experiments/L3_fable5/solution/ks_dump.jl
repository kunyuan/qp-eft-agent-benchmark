# ks_dump.jl — run the pinned KS-DFT setup and dump eigenvalues, Fermi level,
# and plane-wave coefficients of the Bloch states at the explicit k-points
# k_frac = t * endpoint_frac for every t in grid.csv.
#
# Usage:
#   julia --project=<env> ks_dump.jl <element_config.json> <grid.csv> <out_prefix>
#
# Writes:
#   <out_prefix>_meta.json : element, eF_Ha, recip lattice (cols, Bohr^-1),
#                            unit cell volume, atom positions (frac), kpoints,
#                            eigenvalues (Ha) per kpoint, band counts
#   <out_prefix>_psi.bin   : per kpoint: Int64 nG, Int64 nb, Int64 G[3*nG],
#                            ComplexF64 psi[nG*nb] (column-major, band-major last)

using DFTK
using PseudoPotentialData
using JSON
using LinearAlgebra

function main()
    config_path, grid_path, out_prefix = ARGS[1], ARGS[2], ARGS[3]

    cfg = JSON.parsefile(config_path)
    element = cfg["element"]
    st = cfg["structure"]
    dft = cfg["dft"]

    # lattice vectors a1,a2,a3 are given as the three entries of
    # structure.lattice_vectors_bohr; they form the COLUMNS of the lattice matrix
    lat_rows = st["lattice_vectors_bohr"]
    lattice = zeros(Float64, 3, 3)
    for (j, a) in enumerate(lat_rows)   # a = a_j vector
        for i in 1:3
            lattice[i, j] = a[i]
        end
    end

    positions_raw = st["atom_positions_frac"]
    positions = [DFTK.Vec3(Float64.(p)...) for p in positions_raw]
    natoms = length(positions)

    family = PseudoFamily(dft["pseudopotential_family"])
    psp = load_psp(family, Symbol(element))
    atom = ElementPsp(Symbol(element), psp)
    atoms = fill(atom, natoms)

    temperature = Float64(dft["smearing_Ha"])
    model = model_DFT(lattice, atoms, positions;
                      functionals=LDA(),
                      temperature,
                      smearing=DFTK.Smearing.FermiDirac())

    Ecut = Float64(dft["ecut_Ha"])
    kgrid = MonkhorstPack(Int.(dft["kgrid"]))
    basis = PlaneWaveBasis(model; Ecut, kgrid)

    scfres = self_consistent_field(basis; tol=1e-10)
    eF = scfres.εF

    # explicit k-points from the grid
    endpoint = Float64.(st["path"]["endpoint_frac"])
    ts = Float64[]
    point_ids = Int[]
    for (i, line) in enumerate(eachline(grid_path))
        i == 1 && continue
        isempty(strip(line)) && continue
        parts = split(strip(line), ',')
        push!(point_ids, parse(Int, parts[1]))
        push!(ts, parse(Float64, parts[2]))
    end
    kcoords = [DFTK.Vec3((t .* endpoint)...) for t in ts]

    zval = Float64(dft["z_valence"])
    n_electrons = zval * natoms
    n_bands = ceil(Int, n_electrons / 2) + 6

    bands = compute_bands(scfres, ExplicitKpoints(kcoords); n_bands, tol=1e-10)

    bs_basis = bands.basis
    recip = collect(bs_basis.model.recip_lattice)  # columns are b1,b2,b3 in Bohr^-1

    open(out_prefix * "_psi.bin", "w") do io
        for ik in 1:length(kcoords)
            kpt = bs_basis.kpoints[ik]
            Gs = G_vectors(bs_basis, kpt)   # vector of Vec3{Int}
            nG = length(Gs)
            psik = bands.ψ[ik]              # nG x nbands complex
            nb = size(psik, 2)
            write(io, Int64(nG))
            write(io, Int64(nb))
            Gmat = Array{Int64}(undef, 3, nG)
            for (ig, G) in enumerate(Gs)
                Gmat[:, ig] .= Int64.(G)
            end
            write(io, Gmat)
            write(io, ComplexF64.(psik))
        end
    end

    meta = Dict(
        "element" => element,
        "eF_Ha" => eF,
        "recip_lattice_cols_invbohr" => [recip[:, j] for j in 1:3],
        "lattice_cols_bohr" => [lattice[:, j] for j in 1:3],
        "unit_cell_volume_bohr3" => bs_basis.model.unit_cell_volume,
        "atom_positions_frac" => positions_raw,
        "point_ids" => point_ids,
        "ts" => ts,
        "kcoords_frac" => [collect(k) for k in kcoords],
        "eigenvalues_Ha" => [collect(bands.eigenvalues[ik]) for ik in 1:length(kcoords)],
        "n_bands" => n_bands,
        "n_electrons" => n_electrons,
        "smearing_Ha" => temperature,
    )
    open(out_prefix * "_meta.json", "w") do io
        JSON.print(io, meta)
    end
    println("DONE eF_Ha=", eF)
end

main()
