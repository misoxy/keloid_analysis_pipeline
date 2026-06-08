"""Stage 6: all-pair spatial neighbourhood enrichment with spatial-bin null.

Reads:  outputs/<roi_id>/adata_with_igg.h5ad
Writes: outputs/<roi_id>/tables/06_neighbourhood_long.csv
        outputs/<roi_id>/tables/06_significant_pairs_K10.csv
        outputs/<roi_id>/tables/06_IgG_focal_row_K10.csv
        outputs/<roi_id>/figures/06_neighbourhood_K10_heatmap.png

Test:
  For each ordered pair (focal, neighbour), observed neighbour fraction is
  compared to a null built by permuting labels within spatial bins (preserves
  coarse tissue geography). Run at K = 5, 10, 20. BH-FDR across the full matrix.

Run:
    python pipeline/06_neighbourhood.py --config configs/strip_01.yaml
"""
from pathlib import Path
import argparse
import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad

import sys; sys.path.insert(0, str(Path(__file__).parent))
from io_utils import load_config, out_paths, patch_log1p_base
from stats_utils import neighbourhood_enrichment, matrix_to_long
from plot_utils import neighbour_heatmap

PLASMA_LABEL = "IgG_rich_candidate"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = out_paths(cfg)
    n_cfg = cfg["neighbourhood"]
    SEED = 0

    print(f"[stage 6] roi_id   : {cfg['roi_id']}")
    print(f"[stage 6] K list   : {n_cfg['k_list']}")
    print(f"[stage 6] n_perm   : {n_cfg['n_perm']}")
    print(f"[stage 6] bin size : {n_cfg['spatial_bin_size_px']} px (spatial-bin null)")

    a = ad.read_h5ad(paths["with_igg"])
    patch_log1p_base(a)   # defensive: scanpy 1.9 h5ad round-trip workaround
    print(f"[stage 6] loaded   : {a.shape}")

    label_key = "celltype_detailed_v1"
    counts = a.obs[label_key].astype(str).value_counts()
    labels = a.obs[label_key].astype(str).copy()
    rare = counts[counts < n_cfg["min_cells_per_label"]].index
    if len(rare) > 0:
        labels.loc[labels.isin(rare)] = "rare_other"
        print(f"[stage 6] collapsed {len(rare)} rare labels into 'rare_other': {list(rare)}")
    a.obs["celltype_for_neighbourhood"] = pd.Categorical(labels)
    print(f"[stage 6] label set used:\n{a.obs['celltype_for_neighbourhood'].value_counts().to_string()}")

    xy = np.asarray(a.obsm["spatial"])

    all_long = []
    for K in n_cfg["k_list"]:
        t0 = time.time()
        print(f"\n[stage 6] running K={K}, {n_cfg['n_perm']} perms...")
        res = neighbourhood_enrichment(
            xy, a.obs["celltype_for_neighbourhood"].values,
            K=K, n_perm=n_cfg["n_perm"],
            permutation="spatial_bin",
            bin_size_px=n_cfg["spatial_bin_size_px"],
            seed=SEED,
        )
        all_long.append(matrix_to_long(res, K))
        print(f"[stage 6] K={K} done in {time.time()-t0:.1f}s")

        if K == 10:
            neighbour_heatmap(
                res["log2_ratio"], res["qval"],
                title=f"{cfg['roi_id']}: all-pair log2 enrichment K=10 (* q<{n_cfg['alpha']})",
                filename=paths["fig"] / "06_neighbourhood_K10_heatmap.png",
                sig_q=n_cfg["alpha"],
            )

            if PLASMA_LABEL in res["types"]:
                row = pd.DataFrame({
                    "observed": res["observed_fraction"].loc[PLASMA_LABEL],
                    "expected": res["expected_fraction"].loc[PLASMA_LABEL],
                    "log2":     res["log2_ratio"].loc[PLASMA_LABEL],
                    "zscore":   res["zscore"].loc[PLASMA_LABEL],
                    "p_2sided": res["pval"].loc[PLASMA_LABEL],
                    "q_BH_FDR": res["qval"].loc[PLASMA_LABEL],
                }).sort_values("zscore", ascending=False)
                row.to_csv(paths["tab"] / "06_IgG_focal_row_K10.csv")
                print(f"[stage 6] IgG focal row K=10 (top 8):\n{row.head(8).round(3).to_string()}")

    long_df = pd.concat(all_long, ignore_index=True)
    long_df.to_csv(paths["tab"] / "06_neighbourhood_long.csv", index=False)
    print(f"[stage 6] wrote: {paths['tab'] / '06_neighbourhood_long.csv'}")

    sig = long_df[(long_df["K"] == 10) & (long_df["qval"] < n_cfg["alpha"])].copy()
    sig["direction"] = np.where(sig["log2_ratio"] > 0, "enriched", "depleted")
    sig = sig.sort_values(["direction", "log2_ratio"], ascending=[True, False])
    sig.to_csv(paths["tab"] / "06_significant_pairs_K10.csv", index=False)
    print(f"[stage 6] significant pairs K=10 at q<{n_cfg['alpha']}: {len(sig)}")


if __name__ == "__main__":
    main()
