"""Statistical helpers: BH-FDR, marker scoring, neighbourhood enrichment with
spatial-bin permutation, distance label-shuffle test.

All functions are pure: they take arrays in and return arrays / DataFrames out,
no global state. Re-used across stages 4, 6, 7, 8.
"""
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.neighbors import NearestNeighbors


# ============================================================
# FDR
# ============================================================
def benjamini_hochberg(pvals):
    """BH-FDR on any-shape numpy array of p-values. Returns q-values same shape."""
    p = np.asarray(pvals); shape = p.shape
    flat = p.flatten(); n = len(flat)
    order = np.argsort(flat); ranks = np.arange(1, n + 1)
    sorted_p = flat[order]
    raw_adj = sorted_p * n / ranks
    adj_sorted = np.minimum.accumulate(raw_adj[::-1])[::-1]
    adj_sorted = np.minimum(adj_sorted, 1.0)
    adj = np.empty(n); adj[order] = adj_sorted
    return adj.reshape(shape)


# ============================================================
# Gene expression helpers
# ============================================================
def to_1d(x):
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


def gene_vector(adata, gene):
    if gene not in adata.var_names:
        return None
    return to_1d(adata[:, gene].X)


def present_genes(adata, genes):
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    return present, missing


# ============================================================
# Marker scoring (z-scored mean expression of each panel)
# ============================================================
def score_marker_sets(adata, marker_sets, prefix):
    """For each panel: store z-scored mean expression in adata.obs[f'{prefix}_{name}'].

    Returns a presence table so the lab member can see which markers were missing.
    """
    rows = []
    for name, genes in marker_sets.items():
        present, missing = present_genes(adata, genes)
        col = f"{prefix}_{name}"
        if len(present) == 0:
            adata.obs[col] = np.nan
            rows.append({
                "marker_set": name, "n_present": 0, "n_missing": len(missing),
                "present": "", "missing": ",".join(missing),
            })
            continue
        mat = np.vstack([gene_vector(adata, g) for g in present])
        raw_score = mat.mean(axis=0)
        z = (raw_score - raw_score.mean()) / (raw_score.std() + 1e-9)
        adata.obs[col] = z
        rows.append({
            "marker_set": name, "n_present": len(present), "n_missing": len(missing),
            "present": ",".join(present), "missing": ",".join(missing),
        })
    return pd.DataFrame(rows).sort_values(
        ["n_present", "marker_set"], ascending=[False, True])


def marker_fraction_table(adata, label_key, markers):
    """Per-label fraction of cells expressing each marker. Diagnostic table."""
    labels = adata.obs[label_key].astype(str)
    rows = []
    for lab in sorted(labels.unique()):
        mask = (labels == lab).values
        row = {"label": lab, "n_cells": int(mask.sum())}
        for g in markers:
            v = gene_vector(adata, g)
            row[g] = np.nan if v is None else float((v[mask] > 0).mean())
        rows.append(row)
    return pd.DataFrame(rows).set_index("label")


# ============================================================
# Spatial-bin permutation (the conservative null)
# ============================================================
def spatial_bin_permutation(labels, positions, bin_size_px, rng):
    """Permute labels within square spatial bins of side bin_size_px.

    Preserves coarse tissue geography: cells in different parts of the
    tissue are never swapped with each other, only within local bins.

    Returns: permuted labels (same shape as input).
    """
    bin_x = (positions[:, 0] // bin_size_px).astype(int)
    bin_y = (positions[:, 1] // bin_size_px).astype(int)
    bin_id = bin_x * (bin_y.max() + 1) + bin_y
    perm = labels.copy()
    for b in np.unique(bin_id):
        m = bin_id == b
        idx = np.where(m)[0]
        if len(idx) > 1:
            perm[idx] = rng.permutation(labels[idx])
    return perm


# ============================================================
# Neighbourhood enrichment with spatial-bin null
# ============================================================
def neighbourhood_enrichment(positions, labels, K=10, n_perm=500,
                             permutation="spatial_bin", bin_size_px=200, seed=0):
    """Permutation-based pairwise neighbour enrichment.

    permutation: 'spatial_bin' (default, conservative) or 'global' (liberal).

    Returns dict with observed/expected/log2/zscore/pval/qval DataFrames
    indexed by [focal x neighbour] cell types.
    """
    rng = np.random.default_rng(seed)
    labels = np.asarray(labels).astype(str)
    nn = NearestNeighbors(n_neighbors=K + 1).fit(positions)
    _, idx = nn.kneighbors(positions); idx = idx[:, 1:]
    types = np.array(sorted(pd.unique(labels))); nT = len(types)

    def composition(lbls):
        nb = lbls[idx]
        comp = np.zeros((nT, nT))
        for i, a_t in enumerate(types):
            focal = lbls == a_t
            if focal.sum() == 0: continue
            for j, b_t in enumerate(types):
                comp[i, j] = (nb[focal] == b_t).mean()
        return comp

    observed = composition(labels)
    expected = np.tile(
        pd.Series(labels).value_counts(normalize=True).reindex(types).fillna(0).values,
        (nT, 1),
    )

    null = np.zeros((n_perm, nT, nT))
    for p in range(n_perm):
        if permutation == "spatial_bin":
            perm = spatial_bin_permutation(labels, positions, bin_size_px, rng)
        else:
            perm = rng.permutation(labels)
        null[p] = composition(perm)

    zscore = (observed - null.mean(axis=0)) / (null.std(axis=0) + 1e-9)
    p_hi = ((null >= observed[None]).sum(axis=0) + 1) / (n_perm + 1)
    p_lo = ((null <= observed[None]).sum(axis=0) + 1) / (n_perm + 1)
    pval = np.minimum(2 * np.minimum(p_hi, p_lo), 1.0)
    qval = benjamini_hochberg(pval)
    log2_ratio = np.log2((observed + 1e-3) / (expected + 1e-3))

    return {
        "types": types,
        "observed_fraction": pd.DataFrame(observed, index=types, columns=types),
        "expected_fraction": pd.DataFrame(expected, index=types, columns=types),
        "log2_ratio":        pd.DataFrame(log2_ratio, index=types, columns=types),
        "zscore":            pd.DataFrame(zscore, index=types, columns=types),
        "pval":              pd.DataFrame(pval, index=types, columns=types),
        "qval":              pd.DataFrame(qval, index=types, columns=types),
    }


def matrix_to_long(res, K):
    """Tidy a neighbourhood-result dict into long format for CSV export."""
    rows = []
    for f in res["types"]:
        for n in res["types"]:
            rows.append({
                "K": K, "focal": f, "neighbour": n,
                "observed_fraction": res["observed_fraction"].loc[f, n],
                "expected_fraction": res["expected_fraction"].loc[f, n],
                "log2_ratio":        res["log2_ratio"].loc[f, n],
                "zscore":            res["zscore"].loc[f, n],
                "pval":              res["pval"].loc[f, n],
                "qval":              res["qval"].loc[f, n],
            })
    return pd.DataFrame(rows)


# ============================================================
# Distance label-shuffle test
# ============================================================
def distance_label_test(adata, label_key, dist_key, n_perm=500,
                        min_cells=20, permutation="spatial_bin",
                        bin_size_px=200, seed=0):
    """For each label with >= min_cells, test whether observed median distance
    differs from a null (spatial-bin or global shuffle of labels).

    Returns DataFrame: label x [observed/null/delta/pval/qval].
    """
    rng = np.random.default_rng(seed)
    labels = adata.obs[label_key].astype(str).values
    dist = adata.obs[dist_key].values
    positions = np.asarray(adata.obsm["spatial"])
    counts = pd.Series(labels).value_counts()
    types = sorted(counts[counts >= min_cells].index)

    rows = []
    for lab in types:
        mask = labels == lab
        obs_median = float(np.nanmedian(dist[mask]))
        obs_mean = float(np.nanmean(dist[mask]))
        null_median = np.zeros(n_perm); null_mean = np.zeros(n_perm)
        for p in range(n_perm):
            if permutation == "spatial_bin":
                perm = spatial_bin_permutation(labels, positions, bin_size_px, rng)
            else:
                perm = rng.permutation(labels)
            pm = perm == lab
            null_median[p] = np.nanmedian(dist[pm])
            null_mean[p]   = np.nanmean(dist[pm])
        p_hi = ((null_median >= obs_median).sum() + 1) / (n_perm + 1)
        p_lo = ((null_median <= obs_median).sum() + 1) / (n_perm + 1)
        pval = float(min(2 * min(p_hi, p_lo), 1.0))
        rows.append({
            "label": lab, "n_cells": int(mask.sum()),
            "observed_median_dist": obs_median,
            "null_median_dist": float(null_median.mean()),
            "median_delta_vs_null": obs_median - float(null_median.mean()),
            "observed_mean_dist": obs_mean,
            "null_mean_dist": float(null_mean.mean()),
            "pval": pval,
        })
    out = pd.DataFrame(rows)
    out["qval"] = benjamini_hochberg(out["pval"].values.reshape(1, -1)).ravel()
    return out.sort_values("median_delta_vs_null")
