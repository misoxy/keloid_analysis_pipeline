"""Stage 8: local stromal niche around IgG-rich cells (near vs far).

Reads:  outputs/<roi_id>/adata_with_igg.h5ad
Writes: outputs/<roi_id>/tables/08_niche_program_test.csv
        outputs/<roi_id>/tables/08_niche_gene_test.csv
        outputs/<roi_id>/figures/08_near_vs_far_spatial.png

Logic:
  For each stromal cell (fibroblast / deep_dermal_stromal / fibroblast_uncertain),
  compute distance to the nearest IgG-rich candidate cell. Tag as "near" if
  distance < near_threshold_px, "far" if > far_threshold_px. Skip the in-between
  to make the comparison cleaner.

  Then ask: do near-IgG stromal cells score higher on keloid programs (matrix,
  myofibroblast, inflammatory) than far-from-IgG stromal cells? Spatial-bin
  permutation controls for tissue geography.

  This was the analysis that saved the June8 story when the discrete neighbour
  test came back negative - it revealed the IgG aggregate is NOT embedded in an
  activated stromal niche.

Run:
    python pipeline/08_local_niche.py --config configs/strip_01.yaml
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
from stats_utils import (
    benjamini_hochberg, gene_vector, spatial_bin_permutation,
)
from plot_utils import spatial_scatter


STROMAL_LABELS = ["fibroblast", "deep_dermal_stromal", "fibroblast_uncertain"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = out_paths(cfg)
    n_cfg = cfg["local_niche"]
    nb_bin = cfg["neighbourhood"]["spatial_bin_size_px"]
    # Generic: read the focal-population label from the focal block (default to IgG)
    fcfg = cfg.get("focal") or cfg.get("igg") or {}
    FOCAL_LABEL = fcfg.get("target_label", "IgG_rich_candidate")
    SEED = 0

    print(f"[stage 8] roi_id         : {cfg['roi_id']}")
    print(f"[stage 8] focal label    : {FOCAL_LABEL}")
    print(f"[stage 8] near threshold : {n_cfg['near_threshold_px']} px")
    print(f"[stage 8] far threshold  : {n_cfg['far_threshold_px']} px")

    a = ad.read_h5ad(paths["with_igg"])
    patch_log1p_base(a)   # defensive: scanpy 1.9 h5ad round-trip workaround
    xy = np.asarray(a.obsm["spatial"])

    label_key = "celltype_detailed_v1"
    labels = a.obs[label_key].astype(str).values
    is_focal = labels == FOCAL_LABEL
    is_stromal = np.isin(labels, STROMAL_LABELS)
    print(f"[stage 8] focal cells     : {int(is_focal.sum())}")
    print(f"[stage 8] stromal cells   : {int(is_stromal.sum())}")
    if is_focal.sum() == 0:
        print(f"[stage 8] WARNING: no '{FOCAL_LABEL}' cells found. Stage 5 may have been "
              f"disabled or this sample has no focal aggregate. Skipping niche analysis.")
        return

    nn = NearestNeighbors(n_neighbors=1).fit(xy[is_focal])
    d, _ = nn.kneighbors(xy); d = d.ravel()
    a.obs["dist_to_nearest_focal_px"] = d

    near = is_stromal & (d <= n_cfg["near_threshold_px"])
    far  = is_stromal & (d >= n_cfg["far_threshold_px"])
    print(f"[stage 8] near-focal stromal: {int(near.sum())}  "
          f"far-from-focal stromal: {int(far.sum())}")

    if near.sum() < 20 or far.sum() < 20:
        print("[stage 8] WARNING: too few cells in near or far group. Skipping.")
        return

    # ---- Program-level test ----------------------------------------------
    prog_rows = []
    for prog_name, gene_list in n_cfg["niche_programs"].items():
        present = [g for g in gene_list if g in a.var_names]
        if len(present) == 0:
            prog_rows.append({"program": prog_name, "n_present": 0,
                              "near_mean": np.nan, "far_mean": np.nan,
                              "delta": np.nan, "pval": np.nan, "missing": ",".join(gene_list)})
            continue
        mat = np.vstack([gene_vector(a, g) for g in present])
        score = mat.mean(axis=0)
        z = (score - score.mean()) / (score.std() + 1e-9)
        near_mean = float(z[near].mean()); far_mean = float(z[far].mean())
        obs_delta = near_mean - far_mean

        # spatial-bin permutation null: shuffle labels within bins, recompute delta
        null_delta = np.zeros(n_cfg["n_perm"])
        rng = np.random.default_rng(SEED)
        nearfar_label = np.full(len(z), "other", dtype=object)
        nearfar_label[near] = "near"
        nearfar_label[far]  = "far"
        for p in range(n_cfg["n_perm"]):
            perm = spatial_bin_permutation(nearfar_label, xy, nb_bin, rng)
            null_delta[p] = float(z[perm == "near"].mean() - z[perm == "far"].mean())
        p_hi = ((null_delta >= obs_delta).sum() + 1) / (n_cfg["n_perm"] + 1)
        p_lo = ((null_delta <= obs_delta).sum() + 1) / (n_cfg["n_perm"] + 1)
        pval = float(min(2 * min(p_hi, p_lo), 1.0))
        prog_rows.append({
            "program": prog_name, "n_present": len(present),
            "near_mean": near_mean, "far_mean": far_mean,
            "delta": obs_delta, "pval": pval,
            "missing": ",".join([g for g in gene_list if g not in present]),
        })

    prog_df = pd.DataFrame(prog_rows)
    prog_df["qval"] = benjamini_hochberg(prog_df["pval"].fillna(1.0).values.reshape(1, -1)).ravel()
    prog_df = prog_df.sort_values("delta")
    prog_df.to_csv(paths["tab"] / "08_niche_program_test.csv", index=False)
    print(f"[stage 8] niche program test:\n{prog_df.round(3).to_string(index=False)}")

    # ---- Spatial diagnostic plot -----------------------------------------
    nf_label = pd.Series("other", index=a.obs_names, dtype=object)
    nf_label[near] = "stromal_near_focal"
    nf_label[far]  = "stromal_far_focal"
    nf_label[is_focal] = FOCAL_LABEL
    a.obs["near_vs_far_focal"] = pd.Categorical(nf_label)
    spatial_scatter(
        a, "near_vs_far_focal",
        title=f"{cfg['roi_id']}: stromal cells near vs far from {FOCAL_LABEL}",
        filename=paths["fig"] / "08_near_vs_far_spatial.png",
        highlight=[FOCAL_LABEL],
    )

    a.write_h5ad(paths["atlas"])
    print(f"[stage 8] wrote final atlas: {paths['atlas']}")


if __name__ == "__main__":
    main()
