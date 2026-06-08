"""Stage 7: per-cell distance to nearest vessel cell + label-shuffle test.

Reads:  outputs/<roi_id>/adata_with_igg.h5ad
Writes: outputs/<roi_id>/tables/07_distance_test.csv
        outputs/<roi_id>/figures/07_distance_boxplot.png
        outputs/<roi_id>/adata_with_igg.h5ad (in place: add dist_to_nearest_vessel_px)

Test:
  Vessel cells = endothelial + lymphatic_endothelial (configurable).
  For every non-vessel cell, distance in chip pixels to nearest vessel cell.
  Spatial-bin label shuffle: does this label sit closer / farther than chance?

Run:
    python pipeline/07_distance_to_vessel.py --config configs/strip_01.yaml
"""
from pathlib import Path
import argparse
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
from sklearn.neighbors import NearestNeighbors

import sys; sys.path.insert(0, str(Path(__file__).parent))
from io_utils import load_config, out_paths, patch_log1p_base
from stats_utils import distance_label_test
from plot_utils import distance_boxplot


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = out_paths(cfg)
    d_cfg = cfg["distance"]
    nb_bin = cfg["neighbourhood"]["spatial_bin_size_px"]
    SEED = 0

    print(f"[stage 7] roi_id        : {cfg['roi_id']}")
    print(f"[stage 7] vessel labels : {d_cfg['vessel_labels']}")
    print(f"[stage 7] perms         : {d_cfg['n_perm']} (spatial-bin null, bin={nb_bin}px)")

    a = ad.read_h5ad(paths["with_igg"])
    patch_log1p_base(a)   # defensive: scanpy 1.9 h5ad round-trip workaround
    xy = np.asarray(a.obsm["spatial"])
    label_key = "celltype_detailed_v1"

    # Detect vessel cells. Labels in the config might be the broad ("endothelial")
    # or candidate ("blood_endothelial_candidate") variants; check both.
    vessel_label_candidates = set()
    for base in d_cfg["vessel_labels"]:
        vessel_label_candidates.update([base, base + "_candidate", "blood_" + base + "_candidate"])
    actual_vessel_labels = [l for l in vessel_label_candidates if l in a.obs[label_key].astype(str).unique()]
    is_vessel = a.obs[label_key].astype(str).isin(actual_vessel_labels).values
    n_vessel = int(is_vessel.sum())
    print(f"[stage 7] vessel cells found: {n_vessel} (labels: {actual_vessel_labels})")

    if n_vessel == 0:
        print("[stage 7] WARNING: no vessel cells. Skipping distance test.")
        return

    # Distance to nearest vessel
    nn = NearestNeighbors(n_neighbors=1).fit(xy[is_vessel])
    d, _ = nn.kneighbors(xy)
    d = d.ravel()
    d[is_vessel] = np.nan
    a.obs["dist_to_nearest_vessel_px"] = d
    a.obs["dist_to_nearest_vessel_um"] = d / d_cfg["px_per_um"]

    non_v_d = d[~is_vessel]
    print(f"[stage 7] non-vessel dist (px): median={np.nanmedian(non_v_d):.0f}, "
          f"p25-p75={np.nanpercentile(non_v_d,25):.0f}-{np.nanpercentile(non_v_d,75):.0f}, "
          f"max={np.nanmax(non_v_d):.0f}")

    # Label-shuffle test (spatial-bin null)
    test_df = distance_label_test(
        a, label_key=label_key, dist_key="dist_to_nearest_vessel_px",
        n_perm=d_cfg["n_perm"], min_cells=d_cfg["min_cells_per_label"],
        permutation="spatial_bin", bin_size_px=nb_bin, seed=SEED,
    )
    test_df["observed_median_um"] = test_df["observed_median_dist"] / d_cfg["px_per_um"]
    test_df["null_median_um"]     = test_df["null_median_dist"]     / d_cfg["px_per_um"]
    test_df.to_csv(paths["tab"] / "07_distance_test.csv", index=False)
    print(f"[stage 7] label-shuffle test:\n{test_df.round(2).to_string(index=False)}")

    distance_boxplot(
        d, a.obs[label_key].astype(str).values,
        title=f"{cfg['roi_id']}: distance to nearest vessel cell (lower = closer)",
        filename=paths["fig"] / "07_distance_boxplot.png",
        vline=float(np.nanmedian(non_v_d)),
    )

    a.write_h5ad(paths["with_igg"])
    print(f"[stage 7] updated: {paths['with_igg']} (added dist_to_nearest_vessel_px / _um)")


if __name__ == "__main__":
    main()
