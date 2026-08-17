#!/usr/bin/env python3
# Minimal single-cell RNA-seq preprocessing: simulate counts, normalise, and
# select highly variable genes. numpy only (the demo adds sklearn for PCA/kmeans).
import numpy as np


def simulate_counts(n_cells=800, n_genes=1000, n_types=4, markers_per_type=25, seed=0):
    # Each cell belongs to a type; each type over-expresses its own marker block.
    rng = np.random.default_rng(seed)
    labels = rng.integers(0, n_types, n_cells)
    base = rng.lognormal(0.5, 1.0, n_genes)                     # baseline expression per gene
    lam = np.tile(base, (n_cells, 1))
    for t in range(n_types):
        mk = slice(t * markers_per_type, (t + 1) * markers_per_type)
        lam[labels == t, mk] *= rng.uniform(8, 15)              # marker over-expression
    lib = rng.uniform(0.5, 1.5, n_cells)[:, None]               # per-cell sequencing depth
    counts = rng.poisson(lam * lib)
    counts *= (rng.random(counts.shape) > 0.2)                  # dropout
    return counts, labels


def normalise(counts):
    # Library-size normalise to the median depth, then log1p.
    lib = counts.sum(1, keepdims=True) + 1e-9
    cpm = counts / lib * np.median(lib)
    return np.log1p(cpm)


def highly_variable(logn, n_top=200):
    # Pick the most variable genes (the signal-carrying ones).
    v = logn.var(0)
    return np.argsort(v)[::-1][:n_top]


if __name__ == "__main__":
    c, y = simulate_counts()
    ln = normalise(c); hvg = highly_variable(ln)
    print("counts:", c.shape, "median depth:", int(np.median(c.sum(1))), "HVG:", len(hvg))
