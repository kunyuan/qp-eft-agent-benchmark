# gen_packet_data.jl — generate the per-level data files for the agent packet.
#
# For each element, dumps (into reference/packet_data/<El>/):
#   atomic_core_<ns>.csv   r, u_c(r), V_H_c(r)         [L2/L3: agent computes f_c]
#   fc_table_<ns>.csv      K, f_c(K)                   [L1: f_c given]
#   core_model.json        per-channel J_c, f0, DeltaE_c, + meta
#
# Usage: julia --project=<env> gen_packet_data.jl
using Printf, JSON
include(joinpath(@__DIR__, "gold_runner.jl"))   # reuses ELEMENTS, build_form_factors, make_fc_eval

outroot = joinpath(@__DIR__, "packet_data")
mkpath(outroot)

Kgrid = collect(range(0.0, 20.0, length=1001))

for element in ["Na", "Al", "K", "Mg"]
    el = ELEMENTS[element]
    ffs, r, dr = build_form_factors(el)
    eval_fc = make_fc_eval(ffs, r, dr)
    d = joinpath(outroot, element); mkpath(d)

    channels = []
    for ff in ffs
        # trim radial grid to where the core orbital has support
        rcut_i = findlast(i -> abs(ff.u[i]) > 1e-7, eachindex(ff.u))
        rcut_i = something(rcut_i, length(r))
        rng = 1:rcut_i

        open(joinpath(d, "atomic_core_$(ff.name).csv"), "w") do io
            println(io, "r_bohr,u_c,V_H_c")
            for i in rng
                @printf(io, "%.5f,%.8e,%.8e\n", r[i], ff.u[i], ff.VHc[i])
            end
        end
        open(joinpath(d, "fc_table_$(ff.name).csv"), "w") do io
            println(io, "K_bohr_inv,f_c")
            for K in Kgrid
                @printf(io, "%.4f,%.8e\n", K, eval_fc(ff, K))
            end
        end
        push!(channels, Dict("channel"=>ff.name, "J_c_Ha"=>ff.Jc, "f0_Ha"=>ff.f0,
                             "DeltaE_c_Ha"=>ff.ΔE, "r_max_bohr"=>round(r[rcut_i],digits=3)))
        @printf("[%s] %s: ΔE=%.4f Ha, f0=%.4f, J_c=%.4f, r_cut=%.2f, npts=%d\n",
                element, ff.name, ff.ΔE, ff.f0, ff.Jc, r[rcut_i], length(rng))
    end

    open(joinpath(d, "core_model.json"), "w") do io
        JSON.print(io, Dict(
            "element"=>element, "Z_nuc"=>el.Z,
            "core_s_channels"=>channels,
            "form_factor_grid"=>Dict("K_min"=>0.0, "K_max"=>20.0, "n"=>length(Kgrid)),
            "atomic_solver"=>"radial LDA (Dirac exchange), dr=0.002, r_max=40 Bohr; DeltaE_c from ΔSCF",
        ), 2)
    end
end
println("DONE")
