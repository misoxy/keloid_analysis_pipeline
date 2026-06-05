"""Run BANKSY at a sweep of lambda values on one strip.

Reads outputs/<roi_id>/adata_normalised.h5ad, runs BANKSY for each lambda in
the sweep, clusters each result with Leiden, and saves:
  - outputs/<roi_id>/adata_banksy.h5ad
    (the AnnData with one cluster column per lambda, e.g. banksy_l0.2)
  - outputs/<roi_id>/figures/banksy_lambda_sweep.png
    (a panel: UMAP coloured by clusters, one panel per lambda)

Lambda controls how much spatial neighbour information is mixed in:
  - lambda = 0.0 -> ignore neighbours (same as standard Leiden, baseline)
  - lambda ~ 0.2 -> cell typing (the prof's preferred setting)
  - lambda ~ 0.4 -> balanced mix
  - lambda ~ 0.8 -> tissue domain segmentation (epidermal zone vs dermal zone)

Run:
    python pipeline/03_banksy_cluster.py --config configs/strip_01.yaml
"""

from pathlib import Path
import argparse
import warnings
import yaml

import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix

warnings.filterwarnings("ignore")


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def cluster_one_lambda(a_full, banksy_dict, lam, max_m, leiden_resolution):
    """For a single lambda, extract the BANKSY matrix, do PCA -> Leiden -> UMAP.

    Returns the cluster labels (one per cell) and the 2D UMAP coords.
    """
    # banksy_dict is keyed by nbr_weight_decay then by lambda. The inner value
    # is itself an AnnData whose X is the neighbour-augmented BANKSY matrix.
    decay_key = list(banksy_dict.keys())[0]
    banksy_adata = banksy_dict[decay_key][lam]["adata"]

    # Build a fresh AnnData so we don't pollute the input.
    # rows = cells (must match a_full), columns = BANKSY features.
    a = ad.AnnData(X=np.asarray(banksy_adata.X), obs=a_full.obs.copy())
    a.obsm["spatial"] = a_full.obsm["spatial"].copy()

    # Standard scanpy: PCA -> neighbours -> Leiden -> UMAP.
    sc.pp.pca(a, n_comps=20)
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=20)
    sc.tl.leiden(a, resolution=leiden_resolution, random_state=0)
    sc.tl.umap(a, random_state=0)

    return a.obs["leiden"].astype(str), a.obsm["X_umap"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = Path(cfg["output_dir"])
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    in_path  = out_dir / "adata_normalised.h5ad"
    out_path = out_dir / "adata_banksy.h5ad"

    lambdas    = cfg["banksy"]["lambdas"]
    k_geom     = cfg["banksy"]["k_geom"]
    resolution = cfg["banksy"]["leiden_resolution"]
    max_m = 1   # BANKSY default: include AGF (gradient features)

    print(f"roi_id   : {cfg['roi_id']}")
    print(f"lambdas  : {lambdas}")
    print(f"k_geom   : {k_geom}")
    print(f"leiden r : {resolution}")
    print(f"input    : {in_path}")

    a_full = ad.read_h5ad(in_path)
    print(f"loaded   : {a_full.n_obs:,} cells x {a_full.n_vars:,} genes")

    # Restrict to HVGs to speed up BANKSY (it scales with n_genes).
    if "highly_variable" in a_full.var.columns:
        a_hvg = a_full[:, a_full.var["highly_variable"]].copy()
        print(f"using HVGs: {a_hvg.n_vars:,} genes")
    else:
        a_hvg = a_full

    # Step 1: build spatial neighbour graph and compute neighbour features.
    print("\n[1/3] initialising BANKSY (spatial neighbours + features)...")
    banksy_dict = initialize_banksy(
        a_hvg,
        coord_keys=("x", "y", "spatial"),
        num_neighbours=k_geom,
        nbr_weight_decay="scaled_gaussian",
        max_m=max_m,
        plt_edge_hist=False,
        plt_nbr_weights=False,
        plt_agf_angles=False,
        plt_theta=False,
    )

    # Step 2: stack own + neighbour features into a matrix per lambda.
    print("\n[2/3] building BANKSY matrices for each lambda...")
    banksy_dict, _ = generate_banksy_matrix(
        a_hvg, banksy_dict,
        lambda_list=lambdas,
        max_m=max_m,
        plot_std=False,
        verbose=False,
    )

    # Step 3: cluster on each lambda's matrix.
    print(f"\n[3/3] clustering at each lambda (PCA -> Leiden -> UMAP)...")
    all_results = {}
    for lam in lambdas:
        print(f"  lambda = {lam}")
        labels, umap_xy = cluster_one_lambda(
            a_hvg, banksy_dict, lam, max_m, resolution,
        )
        col = f"banksy_l{lam}"
        a_full.obs[col] = labels.values
        a_full.obsm[f"X_umap_l{lam}"] = umap_xy
        all_results[lam] = (labels, umap_xy)
        print(f"    -> {labels.nunique()} clusters")

    # Save the augmented AnnData.
    a_full.write_h5ad(out_path)
    print(f"\nwrote    : {out_path}  ({out_path.stat().st_size/1e6:.1f} MB)")

    # Comparison figure: UMAP per lambda, side-by-side.
    n = len(lambdas)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]
    for ax, lam in zip(axes, lambdas):
        labels, umap_xy = all_results[lam]
        cats = labels.astype("category")
        colors = plt.cm.tab20(np.linspace(0, 1, max(len(cats.cat.categories), 1)))
        for i, c in enumerate(cats.cat.categories):
            mask = (labels == c).values
            ax.scatter(umap_xy[mask, 0], umap_xy[mask, 1],
                       s=2, color=colors[i], alpha=0.7,
                       edgecolors="none", label=c)
        ax.set_title(f"lambda = {lam}  ({labels.nunique()} clusters)")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")

    fig.suptitle(f"BANKSY lambda sweep -- {cfg['roi_id']}  ({a_full.n_obs:,} cells)",
                 fontsize=14)
    plt.tight_layout()
    fig_path = fig_dir / "banksy_lambda_sweep.png"
    plt.savefig(fig_path, dpi=140, bbox_inches="tight")
    print(f"figure   : {fig_path}")


if __name__ == "__main__":
    main()
