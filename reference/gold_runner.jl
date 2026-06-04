# gold_runner.jl — parametrized gold reference for the QP-EFT benchmark.
#
# Given an element's pinned DFT setup + a grid of fractional path coordinates t,
# computes KS bands at k(t)=t*endpoint, applies the coherent frozen-core QP
# correction, and writes predictions. This is the reference solution: it must
# reproduce the authors' Table I and beat bare KS against ARPES.
#
# Usage: julia --project=<env> gold_runner.jl <element> <grid.csv> <out.csv>
using Printf, LinearAlgebra, DFTK, PseudoPotentialData, DelimitedFiles
include(joinpath(@__DIR__, "gold", "sodium", "atomic_hf.jl"))
const Ha2eV = 27.211386245988

# Per-element pinned setup (from authors' eft-psp run_ks.jl + freq_correction.jl)
const ELEMENTS = Dict(
  "Na" => (Z=11, Zval=1, bravais=:bcc, a=8.107, Ecut=15.0, kgrid=[8,8,8],
           endpoint=[0.0,0.0,0.5], # Γ→N
           core=[(1,0,2),(2,0,2),(2,1,6)], full=[(1,0,2),(2,0,2),(2,1,6),(3,0,1)],
           holes=Dict(1=>[(1,0,1),(2,0,2),(2,1,6)], 2=>[(1,0,2),(2,0,1),(2,1,6)])),
)

function prim_lattice(bravais, a)
    bravais == :bcc && return (a/2)*[-1.0 1 1; 1 -1 1; 1 1 -1]
    bravais == :fcc && return (a/2)*[0.0 1 1; 1 0 1; 1 1 0]
    error("unsupported bravais $bravais")
end

# ---- atomic form factors + ΔSCF core excitation energies ----
struct FormFactor; name::String; u::Vector{Float64}; VHc::Vector{Float64}; Jc::Float64; f0::Float64; ΔE::Float64; end

function build_form_factors(el)
    dr=0.002; rmax=40.0
    Ecore = solve_atom(el.Z, Vector{Tuple{Int,Int,Int}}(el.core); dr=dr, r_max=rmax).E_total
    ΔE = Dict{Int,Float64}()
    for (n,cfg) in el.holes
        ΔE[n] = solve_atom(el.Z, Vector{Tuple{Int,Int,Int}}(cfg); dr=dr, r_max=rmax).E_total - Ecore
    end
    atom = solve_atom(el.Z, Vector{Tuple{Int,Int,Int}}(el.full); dr=dr, r_max=rmax)
    r=atom.rgrid; N=length(r)
    ffs = FormFactor[]
    for orb in atom.orbitals[1:end-1]   # skip valence
        orb.l != 0 && continue
        VHc = orbital_coulomb_potential(orb.u, 1.0, r, dr)
        Jc  = dr*sum(orb.u[i]^2*VHc[i] for i in 1:N)
        Sc  = dr*sum(r[i]*orb.u[i] for i in 1:N)
        Wc  = dr*sum(r[i]*orb.u[i]*VHc[i] for i in 1:N)
        f0  = sqrt(4π)*(Wc - Jc*Sc)
        push!(ffs, FormFactor("$(orb.n)s", copy(orb.u), VHc, Jc, f0, ΔE[orb.n]))
    end
    return ffs, r, dr
end

function make_fc_eval(ffs, r, dr)
    N=length(r)
    function eval_fc(ff::FormFactor, K)
        K < 1e-10 && return ff.f0
        s=0.0; @inbounds for i in 1:N; s += ff.u[i]*(ff.VHc[i]-ff.Jc)*sin(K*r[i]); end
        return sqrt(4π)/K*s*dr
    end
    return eval_fc
end

function run(element, gridfile, outfile)
    el = ELEMENTS[element]
    # 1) atomic side
    ffs, r, dr = build_form_factors(el)
    eval_fc = make_fc_eval(ffs, r, dr)
    @printf("[%s] core channels: %s\n", element,
            join(["$(ff.name) f0=$(round(ff.f0,digits=4)) ΔE=$(round(ff.ΔE,digits=4)) Δc0=$(round(ff.f0^2/ff.ΔE^2,digits=4))" for ff in ffs], " | "))

    # 2) KS side
    lattice = prim_lattice(el.bravais, el.a)
    psp = load_psp(PseudoFamily("cp2k.nc.sr.lda.v0_1.largecore.gth"), Symbol(element))
    model = model_LDA(lattice, [ElementPsp(Symbol(element), psp)], [zeros(3)];
                      temperature=0.001, smearing=DFTK.Smearing.FermiDirac())
    basis = PlaneWaveBasis(model; Ecut=el.Ecut, kgrid=el.kgrid)
    scfres = self_consistent_field(basis; tol=1e-8, mixing=KerkerMixing(),
                                   is_converged=DFTK.ScfConvergenceEnergy(1e-8), callback=identity)
    εF = scfres.εF

    # 3) grid → explicit k(t)
    grid = readdlm(gridfile, ','; header=true)[1]
    pids = Int.(grid[:,1]); ts = Float64.(grid[:,2])
    kpts = [t .* el.endpoint for t in ts]
    bands = compute_bands(scfres, DFTK.ExplicitKpoints(kpts); n_bands=4)
    recip = model.recip_lattice

    # 4) coherent correction, lowest band (occupied)
    rows = Vector{Any}()
    for (i,t) in enumerate(ts)
        kfrac = kpts[i]
        Gs = collect(DFTK.G_vectors(bands.basis, bands.basis.kpoints[i]))
        for n in 1:1   # lowest band (occupied for these simple metals along the path)
            ψ = bands.ψ[i][:,n]
            Δ = 0.0
            for ff in ffs
                Fc = 0.0 + 0.0im
                @inbounds for (ig,G) in enumerate(Gs)
                    K = norm(recip*(kfrac .+ G))
                    Fc += ψ[ig]*eval_fc(ff,K)/abs(ff.ΔE)
                end
                Δ += abs2(Fc)
            end
            eks = (bands.eigenvalues[i][n]-εF)*Ha2eV
            eqp = eks/(1+Δ)   # (E-εF) scaled by z_core=1/(1+Δ)
            push!(rows, (element, pids[i], t, eqp, eks, Δ))
        end
    end

    open(outfile,"w") do io
        println(io,"element,point_id,t,E_pred_eV,E_KS_eV,Delta")
        for r in rows
            @printf(io,"%s,%d,%.6f,%.5f,%.5f,%.6f\n", r...)
        end
    end
    @printf("[%s] εF=%.4f eV, wrote %d rows → %s\n", element, εF*Ha2eV, length(rows), outfile)
end

run(ARGS[1], ARGS[2], ARGS[3])
