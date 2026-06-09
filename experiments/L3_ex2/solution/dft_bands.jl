# dft_bands.jl  -- run one SCF and emit KS bands + plane-wave coefficients at the
# requested k-points.  Called by run_qp.py.  Writes a JSON with, per k-point:
#   eps (eigenvalues, Ha), epsF (Ha), and for each band the |k+G| (Bohr^-1) and
#   the complex coefficients c_nk(G).  Output is intentionally verbose JSON so the
#   Python side does all physics.
#
# Usage: julia dft_bands.jl <config.json> <grid.csv> <out.json> <n_bands>

using DFTK
using PseudoPotentialData
using JSON
using LinearAlgebra

function main()
    config_path = ARGS[1]
    grid_path   = ARGS[2]
    out_path    = ARGS[3]
    n_bands     = parse(Int, ARGS[4])

    cfg = JSON.parsefile(config_path)
    element = cfg["element"]
    st = cfg["structure"]

    # lattice vectors are the rows in JSON; they are a1,a2,a3 -> COLUMNS of lattice
    lv = st["lattice_vectors_bohr"]            # vector of 3 vectors
    lattice = zeros(3,3)
    for j in 1:3, i in 1:3
        lattice[i,j] = lv[j][i]                # column j = a_j
    end

    dft = cfg["dft"]
    family = dft["pseudopotential_family"]
    Ecut = Float64(dft["ecut_Ha"])
    kg = dft["kgrid"]
    kgrid = (Int(kg[1]), Int(kg[2]), Int(kg[3]))
    smearing_Ha = Float64(dft["smearing_Ha"])

    psp = load_psp(PseudoFamily(family), Symbol(element))
    el_atom = ElementPsp(Symbol(element), psp)

    positions = Vector{Vector{Float64}}()
    for p in st["atom_positions_frac"]
        push!(positions, Float64.(p))
    end
    atoms = [el_atom for _ in positions]

    model = model_DFT(lattice, atoms, positions;
                      functionals=[:lda_x, :lda_c_pw],
                      temperature=smearing_Ha,
                      smearing=DFTK.Smearing.FermiDirac())
    basis = PlaneWaveBasis(model; Ecut=Ecut, kgrid=kgrid)
    scfres = self_consistent_field(basis; tol=1e-8)

    # endpoint and ts
    endpoint = Float64.(st["path"]["endpoint_frac"])
    # read grid.csv: header point_id,t
    ts = Float64[]
    pids = Int[]
    open(grid_path) do io
        header = readline(io)
        for line in eachline(io)
            line = strip(line)
            isempty(line) && continue
            parts = split(line, ",")
            push!(pids, parse(Int, strip(parts[1])))
            push!(ts, parse(Float64, strip(parts[2])))
        end
    end

    kpts_frac = [t .* endpoint for t in ts]
    b = compute_bands(scfres, DFTK.ExplicitKpoints(kpts_frac); n_bands=n_bands)

    result = Dict{String,Any}()
    result["epsF"] = b.εF
    result["element"] = element
    kdata = Vector{Any}()
    for ik in 1:length(ts)
        eps = b.eigenvalues[ik]                       # vector length n_bands (Ha)
        kpg = norm.(DFTK.Gplusk_vectors_cart(b.basis, b.basis.kpoints[ik]))  # |k+G| Bohr^-1
        psi = b.ψ[ik]                                 # (n_G x n_bands) complex
        nG = size(psi,1)
        nb = size(psi,2)
        bands = Vector{Any}()
        for n in 1:nb
            cre = real.(psi[:,n])
            cim = imag.(psi[:,n])
            push!(bands, Dict("eps"=>eps[n], "c_re"=>cre, "c_im"=>cim))
        end
        push!(kdata, Dict("point_id"=>pids[ik], "t"=>ts[ik],
                          "kpg"=>kpg, "bands"=>bands, "nG"=>nG, "nb"=>nb))
    end
    result["kpoints"] = kdata

    open(out_path, "w") do io
        JSON.print(io, result)
    end
    println("WROTE ", out_path)
end

main()
