# Pinned environment for the QP-EFT benchmark: Julia + DFTK 0.7.25 + Python.
# No PackageCompiler sysimage — for a single-agent exploration harness the
# one-time precompile/JIT cost is acceptable. DFTK is version-pinned via
# Manifest.toml (gold/thresholds are calibrated to 0.7.25) and the GTH pseudo
# artifact is pre-fetched so the runtime needs no network.
FROM julia:1.12.1-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-numpy python3-scipy ca-certificates && \
    rm -rf /var/lib/apt/lists/*

ENV JULIA_PROJECT=/opt/jenv \
    JULIA_DEPOT_PATH=/opt/julia_depot

WORKDIR /opt/jenv
COPY Project.toml Manifest.toml ./
RUN julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' && \
    julia -e 'using PseudoPotentialData; fam = PseudoFamily("cp2k.nc.sr.lda.v0_1.largecore.gth"); foreach(el -> fam[el], [:Na, :Al, :K, :Mg, :Ca, :Li, :Si])'

WORKDIR /app
COPY packet/ /app/packet/
