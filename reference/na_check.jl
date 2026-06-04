using DFTK, PseudoPotentialData, Unitful, UnitfulAtomic, LinearAlgebra
using Brillouin

a = 8.11  # bohr, bcc conventional lattice constant
lattice = (a/2) * [ -1.0  1.0  1.0;
                     1.0 -1.0  1.0;
                     1.0  1.0 -1.0 ]   # columns = primitive bcc vectors
fam = PseudoFamily("cp2k.nc.sr.lda.v0_1.largecore.gth")
Na = ElementPsp(:Na, fam)
atoms = [Na]; positions = [zeros(3)]

model = model_DFT(lattice, atoms, positions;
                  functionals=LDA(),
                  temperature=0.001, smearing=Smearing.FermiDirac())
basis = PlaneWaveBasis(model; Ecut=30, kgrid=(16,16,16))
scfres = self_consistent_field(basis; tol=1e-8)
εF = scfres.εF
println("εF (Ha) = ", εF)

# bcc N point (Setyawan-Curtarolo) in fractional reciprocal coords
N = [0.0, 0.0, 0.5]
Γ = [0.0, 0.0, 0.0]
ts = range(0, 1; length=21)
kpath = [Γ .+ t .* (N .- Γ) for t in ts]

bands = compute_bands(scfres, DFTK.ExplicitKpoints(kpath); n_bands=4)
ev = bands.eigenvalues  # vector per kpoint
println("\n t      E_band1 - εF (eV)")
for (i,t) in enumerate(ts)
    e1 = (ev[i][1] - εF) * 27.211386245988
    println(rpad(round(t,digits=3),7), round(e1, digits=4))
end
eΓ = (ev[1][1] - εF)*27.211386245988
println("\nKS Γ-point depth (E(Γ)-εF) = ", round(eΓ,digits=3), " eV   [paper: -3.27 eV]")
println("z_core=0.80 → QP Γ depth = ", round(0.80*eΓ,digits=3), " eV   [paper: -2.62 eV, ARPES -2.65..-2.78]")
