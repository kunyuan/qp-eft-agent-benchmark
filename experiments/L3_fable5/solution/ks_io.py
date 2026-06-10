"""Read the ks_dump.jl output (meta json + psi binary)."""
import json
import numpy as np


def read_dump(prefix):
    """Returns (meta, kdata) where kdata[ik] = dict(G=int (nG,3), psi=(nG,nb) complex)."""
    with open(prefix + "_meta.json") as f:
        meta = json.load(f)
    nk = len(meta["kcoords_frac"])
    kdata = []
    with open(prefix + "_psi.bin", "rb") as f:
        for _ in range(nk):
            nG = int(np.fromfile(f, dtype=np.int64, count=1)[0])
            nb = int(np.fromfile(f, dtype=np.int64, count=1)[0])
            # julia wrote Gmat (3, nG) column-major -> consecutive triples are G vectors
            G = np.fromfile(f, dtype=np.int64, count=3 * nG).reshape(nG, 3)
            # psik (nG, nb) column-major -> consecutive nG block per band
            psi = np.fromfile(f, dtype=np.complex128, count=nG * nb).reshape(nb, nG).T
            kdata.append({"G": G, "psi": psi})
    return meta, kdata
