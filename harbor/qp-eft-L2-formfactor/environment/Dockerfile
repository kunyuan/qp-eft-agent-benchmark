# Pinned environment for the QP-EFT benchmark: Julia + DFTK 0.7.25 + Python.
# We do NOT build a PackageCompiler sysimage — for a single-agent exploration
# harness the one-time precompile/JIT cost is acceptable. We DO pin DFTK's
# version via Manifest.toml (the gold/thresholds are calibrated to 0.7.25) and
# pre-fetch the pseudopotential artifact so the runtime needs no network.
FROM julia:1.12.1-bookworm

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-numpy python3-scipy ca-certificates && \
    rm -rf /var/lib/apt/lists/*

ENV JULIA_PROJECT=/opt/jenv \
    JULIA_DEPOT_PATH=/opt/julia_depot

WORKDIR /opt/jenv
COPY Project.toml Manifest.toml ./
# instantiate (installs + precompiles the pinned deps into the baked depot),
# then trigger the GTH pseudo-family artifact download for the benchmark elements
RUN julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' && \
    julia -e 'using PseudoPotentialData; fam = PseudoFamily("cp2k.nc.sr.lda.v0_1.largecore.gth"); foreach(el -> fam[el], [:Na, :Al, :K, :Mg, :Ca, :Li, :Si])'

# Agent working directory + public development packet (agent-facing).
WORKDIR /app
COPY packet/ /app/packet/
