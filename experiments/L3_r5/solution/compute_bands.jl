# compute_bands.jl
# One SCF, bands at explicit k = t*endpoint_frac, and per band a radial profile
# of |psi_nk(r)|^2 around the atoms, obtained from the REAL-SPACE wavefunction on
# DFTK's FFT grid (fast & exact to grid resolution).
#
# We need <psi_nk | F_c | psi_nk> = int_cell |psi_nk(r)|^2 F_c(|r-R_a|) d3r, with
# F_c short-ranged about each atom R_a. We compute it by binning the real-space
# grid points by distance to the NEAREST atom (minimum image), producing a
# histogram  P(r_bin) = sum over grid pts in bin of |psi|^2 * dV.  Then
#   <F_c>_nk = sum_bins P(r_bin) * F_c(r_bin).
# This is exact (no spherical-average approximation): every cell grid point is
# assigned the F_c of its nearest atom. dV = Omega / Nx/Ny/Nz.

using DFTK
using PseudoPotentialData
using JSON
using LinearAlgebra

function main()
    specfile = ARGS[1]
    spec = JSON.parsefile(specfile)

    element = spec["element"]
    A = spec["lattice_vectors_bohr"]
    lattice = zeros(3,3)
    for i in 1:3, j in 1:3
        lattice[j, i] = A[i][j]            # a_i is the i-th COLUMN
    end
    positions = [Float64.(p) for p in spec["atom_positions_frac"]]

    Ecut = Float64(spec["ecut_Ha"])
    kgrid = Int.(spec["kgrid"])
    smearing_Ha = Float64(spec["smearing_Ha"])
    endpoint = Float64.(spec["endpoint_frac"])
    ts = Float64.(spec["ts"])
    n_bands = Int(spec["n_bands"])
    rbins = Float64.(spec["rbins"])        # bin EDGES (Bohr), ascending

    psp = load_psp(PseudoFamily(spec["pseudopotential_family"]), Symbol(element))
    atoms = [ElementPsp(Symbol(element), psp) for _ in positions]

    model = model_DFT(lattice, atoms, positions;
                      functionals=[:lda_x, :lda_c_pw],
                      temperature=smearing_Ha,
                      smearing=DFTK.Smearing.FermiDirac())
    basis = PlaneWaveBasis(model; Ecut=Ecut, kgrid=kgrid)
    scfres = self_consistent_field(basis; tol=1e-8)
    eF = scfres.εF

    kpts = [t .* endpoint for t in ts]
    b = compute_bands(scfres, DFTK.ExplicitKpoints(kpts); n_bands=n_bands)

    Omega = model.unit_cell_volume
    fft_size = b.basis.fft_size
    Ngrid = prod(fft_size)
    dV = Omega / Ngrid

    # Real-space fractional grid points r_frac and their CARTESIAN distance to the
    # nearest atom under minimum image. Precompute the bin index of each grid pt.
    nbin = length(rbins) - 1
    rmax = rbins[end]
    # fractional grid
    Nx, Ny, Nz = fft_size
    # Precompute nearest-atom distance -> bin index for every grid point.
    binidx = Vector{Int}(undef, Ngrid)  # 0 means "outside last bin" (ignored)
    lat = model.lattice
    # neighbor shifts for minimum image
    shifts = [(i,j,k) for i in -1:1 for j in -1:1 for k in -1:1]
    @inbounds for iz in 0:Nz-1, iy in 0:Ny-1, ix in 0:Nx-1
        idx = ix + Nx*(iy + Ny*iz) + 1   # column-major, ix fastest (matches vec(array))
        f = (ix/Nx, iy/Ny, iz/Nz)
        dmin = Inf
        for R in positions
            df0 = (f[1]-R[1], f[2]-R[2], f[3]-R[3])
            for s in shifts
                df = (df0[1]+s[1], df0[2]+s[2], df0[3]+s[3])
                # cartesian
                cx = lat[1,1]*df[1]+lat[1,2]*df[2]+lat[1,3]*df[3]
                cy = lat[2,1]*df[1]+lat[2,2]*df[2]+lat[2,3]*df[3]
                cz = lat[3,1]*df[1]+lat[3,2]*df[2]+lat[3,3]*df[3]
                d = sqrt(cx*cx+cy*cy+cz*cz)
                if d < dmin; dmin = d; end
            end
        end
        if dmin >= rmax
            binidx[idx] = 0
        else
            # find bin: rbins ascending
            bi = searchsortedlast(rbins, dmin)
            binidx[idx] = (bi >= 1 && bi <= nbin) ? bi : 0
        end
    end

    results = Dict{String,Any}()
    results["eF_Ha"] = eF
    results["Omega"] = Omega
    results["dV"] = dV
    kp_out = Vector{Any}()

    for (ik, t) in enumerate(ts)
        eigs = b.eigenvalues[ik]
        psi = b.ψ[ik]
        kpt = b.basis.kpoints[ik]
        nb = size(psi, 2)
        prof = zeros(nb, nbin)
        for bnd in 1:nb
            psir = DFTK.ifft(b.basis, kpt, psi[:, bnd])  # real-space, normalized over cell
            absr = abs2.(psir)                            # |psi|^2
            v = vec(absr)
            @inbounds for gi in 1:Ngrid
                bi = binidx[gi]
                if bi != 0
                    prof[bnd, bi] += v[gi]
                end
            end
        end
        prof .*= dV   # now prof[bnd, bi] = int_{bin} |psi|^2 dV = probability mass in shell bin
        push!(kp_out, Dict(
            "t" => t,
            "eigenvalues_Ha" => collect(eigs),
            "prof" => [collect(prof[bnd, :]) for bnd in 1:nb],
        ))
    end
    results["kpoints"] = kp_out
    results["rbins"] = rbins

    open(spec["out_json"], "w") do io
        JSON.print(io, results)
    end
    println("WROTE ", spec["out_json"])
end

main()
