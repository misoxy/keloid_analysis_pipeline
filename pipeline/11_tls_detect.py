"""Stage 11: tertiary lymphoid structure (TLS) detection.

Reads:  outputs/<roi_id>/adata_with_igg.h5ad
        configs/marker_panels.yaml (tls_markers panel)
Writes: outputs/<roi_id>/tables/11_tls_candidate_cells.csv
        outputs/<roi_id>/tables/11_tls_clusters.csv
        outputs/<roi_id>/figures/11_tls_spatial.png

Logic:
  TLS = ectopic lymphoid aggregate. In skin / keloid Stereo-seq, the candidate
  is a spatial cluster of cells co-expressing B cell, follicular helper T,
  follicular dendritic, and IgG-producing programs. The pipeline:

  1. Score each cell against TLS marker panels (FDC, GC B, Tfh, HEV).
  2. Define "TLS candidate" cells as those with a high TLS score on ANY panel
     OR labelled b_cell_candidate / IgG_rich_candidate / dendritic_cell.
  3. Cluster candidate cells spatially with DBSCAN (eps default 100 px,
     min_samples 10). Each cluster is a candidate TLS.
  4. Report per-cluster: cell count, area, composition by celltype_detailed_v1.

Run:
    python pipeline/11_tls_detect.py --config configs/strip_01.yaml
"""
from pathlib import Path
import argparse
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt

import sys; sys.path.insert(0, str(Path(__file__).parent))
from io_utils import load_config, load_marker_panels, out_paths, patch_log1p_base
from stats_utils import score_marker_sets


TLS_RELEVANT_LABELS = {
    "b_cell_candidate", "IgG_rich_candidate", "dendritic_cell",
    "t_cell_candidate", "lymphatic_endothelial_candidate",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--score-threshold", type=float, default=1.0,
                    help="Z-score above which a cell counts as a TLS candidate")
    ap.add_argument("--dbscan-eps-px", type=float, default=100,
                    help="DBSCAN neighbourhood radius (chip pixels)")
    ap.add_argument("--dbscan-min-samples", type=int, default=10,
                    help="DBSCAN min cluster size")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = out_paths(cfg)
    panels_path = cfg["annotation"]["marker_panels"]
    panels = load_marker_panels(panels_path)

    if "tls_markers" not in panels:
        print("[stage 11] no tls_markers panel in marker_panels.yaml; skipping.")
        return

    print(f"[stage 11] roi_id    : {cfg['roi_id']}")
    print(f"[stage 11] z thresh  : {args.score_threshold}")
    print(f"[stage 11] DBSCAN    : eps={args.dbscan_eps_px}px, min_samples={args.dbscan_min_samples}")

    a = ad.read_h5ad(paths["with_igg"])
    patch_log1p_base(a)
    print(f"[stage 11] loaded   : {a.shape}")

    tls_info = score_marker_sets(a, panels["tls_markers"], prefix="tls")
    print(f"[stage 11] panel coverage:\n{tls_info[['marker_set', 'n_present', 'n_missing']].to_string(index=False)}")

    score_cols = [c for c in [f"tls_{k}" for k in panels["tls_markers"]] if c in a.obs.columns]
    if len(score_cols) == 0:
        print("[stage 11] no TLS panels had any present genes; skipping.")
        return

    max_tls_score = a.obs[score_cols].max(axis=1)
    is_high_score = (max_tls_score >= args.score_threshold).values

    label_col = "celltype_detailed_v1" if "celltype_detailed_v1" in a.obs.columns else None
    is_relevant_label = np.zeros(a.n_obs, dtype=bool)
    if label_col is not None:
        is_relevant_label = a.obs[label_col].astype(str).isin(TLS_RELEVANT_LABELS).values

    is_candidate = is_high_score | is_relevant_label
    n_cand = int(is_candidate.sum())
    print(f"[stage 11] TLS candidate cells: {n_cand} "
          f"({int(is_high_score.sum())} by score, {int(is_relevant_label.sum())} by label)")

    if n_cand < args.dbscan_min_samples:
        print("[stage 11] too few candidates for DBSCAN; nothing to report.")
        return

    xy = np.asarray(a.obsm["spatial"])
    db = DBSCAN(eps=args.dbscan_eps_px, min_samples=args.dbscan_min_samples)
    cluster_id = np.full(a.n_obs, -1, dtype=int)
    cluster_id[is_candidate] = db.fit_predict(xy[is_candidate])
    a.obs["tls_cluster_id"] = cluster_id

    cand_df = pd.DataFrame({
        "cell": a.obs_names,
        "x": xy[:, 0], "y": xy[:, 1],
        "tls_cluster_id": cluster_id,
        "max_tls_score": max_tls_score.values,
        "label": a.obs[label_col].astype(str).values if label_col else "",
        "is_candidate": is_candidate,
    })
    cand_df.to_csv(paths["tab"] / "11_tls_candidate_cells.csv", index=False)

    cluster_rows = []
    for cid in sorted(set(cluster_id)):
        if cid < 0:
            continue
        m = cluster_id == cid
        comp = (a.obs.loc[m, label_col].astype(str).value_counts()
                if label_col else pd.Series(dtype=int))
        cluster_rows.append({
            "tls_cluster_id": cid,
            "n_cells": int(m.sum()),
            "x_centroid": float(xy[m, 0].mean()),
            "y_centroid": float(xy[m, 1].mean()),
            "x_extent_px": float(xy[m, 0].max() - xy[m, 0].min()),
            "y_extent_px": float(xy[m, 1].max() - xy[m, 1].min()),
            "label_composition": "; ".join(f"{lab}={n}" for lab, n in comp.items()),
        })
    cluster_df = pd.DataFrame(cluster_rows).sort_values("n_cells", ascending=False)
    cluster_df.to_csv(paths["tab"] / "11_tls_clusters.csv", index=False)
    print(f"[stage 11] candidate TLS clusters detected: {len(cluster_df)}")
    if len(cluster_df) > 0:
        print(f"[stage 11] top clusters:\n{cluster_df.head(10).to_string(index=False)}")

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.scatter(xy[:, 0], xy[:, 1], s=1, c="lightgrey", alpha=0.25, edgecolors="none")
    ax.scatter(xy[is_candidate & (cluster_id < 0), 0],
               xy[is_candidate & (cluster_id < 0), 1],
               s=6, c="lightblue", alpha=0.6, edgecolors="none",
               label="TLS candidate (noise)")
    palette = plt.cm.tab20(np.linspace(0, 1, max(len(cluster_df), 1)))
    for i, cid in enumerate(cluster_df["tls_cluster_id"]):
        m = cluster_id == cid
        ax.scatter(xy[m, 0], xy[m, 1], s=14, color=palette[i % len(palette)],
                   edgecolors="black", linewidths=0.3, alpha=0.9,
                   label=f"TLS_{int(cid):02d} (n={int(m.sum())})")
    ax.invert_yaxis(); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{cfg['roi_id']}: candidate tertiary lymphoid structures "
                 f"(DBSCAN eps={args.dbscan_eps_px}px, min_samples={args.dbscan_min_samples})")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=8, markerscale=2)
    plt.tight_layout()
    plt.savefig(paths["fig"] / "11_tls_spatial.png", dpi=180, bbox_inches="tight")
    plt.close()
    print(f"[stage 11] wrote: {paths['fig'] / '11_tls_spatial.png'}")


if __name__ == "__main__":
    main()
