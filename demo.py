#!/usr/bin/env python3
# Self-contained demo: a miniature single-cell RNA-seq pipeline, the OSCA
# workflow in one file. Simulate a count matrix with hidden cell types, then
# normalise -> select highly variable genes -> PCA -> embed -> cluster ->
# find marker genes, and check the clusters against the true cell types.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.metrics import adjusted_rand_score
from scrna import simulate_counts, normalise, highly_variable

N_TYPES = 4


def main():
    os.makedirs("figures", exist_ok=True); os.makedirs("results", exist_ok=True)
    counts, truth = simulate_counts(n_cells=800, n_genes=1000, n_types=N_TYPES, seed=0)
    logn = normalise(counts)
    hvg = highly_variable(logn, 200)
    X = StandardScaler().fit_transform(logn[:, hvg])
    pcs = PCA(30, random_state=0).fit(X)
    emb30 = pcs.transform(X)
    emb2 = TSNE(2, perplexity=30, init="pca", random_state=0).fit_transform(emb30)
    clusters = KMeans(N_TYPES, n_init=10, random_state=0).fit_predict(emb30)
    ari = adjusted_rand_score(truth, clusters)

    # Marker genes: top genes by cluster-vs-rest mean difference.
    markers = {}
    for c in range(N_TYPES):
        d = logn[clusters == c].mean(0) - logn[clusters != c].mean(0)
        markers[c] = np.argsort(d)[::-1][:3]

    fig, ax = plt.subplots(2, 2, figsize=(13, 10))

    # Panel 1: embedding coloured by discovered cluster.
    a = ax[0, 0]
    for c in range(N_TYPES):
        m = clusters == c; a.scatter(emb2[m, 0], emb2[m, 1], s=8, label=f"cluster {c}")
    a.set_title(f"t-SNE coloured by cluster (ARI vs truth = {ari:.2f})"); a.set_xticks([]); a.set_yticks([]); a.legend(fontsize=8)

    # Panel 2: embedding coloured by true cell type.
    a = ax[0, 1]
    for t in range(N_TYPES):
        m = truth == t; a.scatter(emb2[m, 0], emb2[m, 1], s=8, label=f"type {t}")
    a.set_title("Same embedding coloured by true cell type"); a.set_xticks([]); a.set_yticks([]); a.legend(fontsize=8)

    # Panel 3: marker-gene mean expression per cluster (heatmap).
    a = ax[1, 0]
    gene_ids = sum([list(markers[c]) for c in range(N_TYPES)], [])
    M = np.array([[logn[clusters == c][:, g].mean() for g in gene_ids] for c in range(N_TYPES)])
    Mz = (M - M.mean(0)) / (M.std(0) + 1e-9)
    im = a.imshow(Mz, aspect="auto", cmap="viridis")
    a.set_xticks(range(len(gene_ids))); a.set_xticklabels([f"g{g}" for g in gene_ids], rotation=90, fontsize=6)
    a.set_yticks(range(N_TYPES)); a.set_yticklabels([f"cluster {c}" for c in range(N_TYPES)])
    a.set_title("Marker genes (z-scored mean expression)"); fig.colorbar(im, ax=a, fraction=0.046)

    # Panel 4: variance explained by the top PCs.
    a = ax[1, 1]
    ve = pcs.explained_variance_ratio_[:20] * 100
    a.bar(range(1, 21), ve, color="#4C72B0")
    a.set_xlabel("principal component"); a.set_ylabel("% variance explained")
    a.set_title("PCA scree: signal in the first few PCs")

    fig.suptitle("Single-cell RNA-seq clustering & marker genes (synthetic data)", fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig("figures/demo.png", dpi=120)

    rows = [{"cluster": c, "marker_genes": ";".join(f"g{g}" for g in markers[c])} for c in range(N_TYPES)]
    pd.DataFrame(rows).to_csv("results/summary.csv", index=False)
    print(f"Cells: {counts.shape[0]}, genes: {counts.shape[1]}, clusters: {N_TYPES}")
    print(f"Adjusted Rand Index (clusters vs true types): {ari:.3f}")
    print("Wrote figures/demo.png and results/summary.csv")


if __name__ == "__main__":
    main()
