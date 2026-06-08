"""Stage 10: aggregate results across multiple ROIs.

Reads:  outputs/<roi_id>/tables/06_IgG_focal_row_K10.csv  for each ROI
        outputs/<roi_id>/adata_atlas_full.h5ad             (optional, for counts)
Writes: outputs/aggregate/cross_roi_IgG_neighbour_zscore.csv
        outputs/aggregate/cross_roi_IgG_neighbour_correlations.csv
        outputs/aggregate/cross_roi_celltype_counts.csv
        outputs/aggregate/cross_roi_summary.html

Logic:
  - For each ROI, pull the IgG focal-row neighbour z-scores at K=10.
  - Build a (ROI x neighbour) matrix of z-scores.
  - Compute pairwise Spearman rank correlation between ROIs. A correlation
    above ~0.5 between two ROIs means the IgG-neighbour pattern reproduces.
  - Mean and stability (sign-consistency) of each neighbour z-score across ROIs.
  - Aggregate cell-type counts per ROI in one table.

Run:
    python pipeline/10_aggregate_rois.py --roi-glob 'outputs/autotile_*' --out outputs/aggregate
"""
from pathlib import Path
import argparse
import glob
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt


def collect_igg_rows(roi_dirs):
    rows = {}
    for roi_dir in roi_dirs:
        roi = Path(roi_dir).name
        csv = Path(roi_dir) / "tables" / "06_IgG_focal_row_K10.csv"
        if not csv.exists():
            continue
        df = pd.read_csv(csv, index_col=0)
        rows[roi] = df["zscore"]
    if not rows:
        return None
    return pd.DataFrame(rows).fillna(0.0)


def collect_celltype_counts(roi_dirs):
    rows = []
    for roi_dir in roi_dirs:
        roi = Path(roi_dir).name
        h5ad = Path(roi_dir) / "adata_atlas_full.h5ad"
        if not h5ad.exists():
            continue
        a = ad.read_h5ad(h5ad)
        col = "celltype_detailed_v1" if "celltype_detailed_v1" in a.obs.columns else None
        if col is None:
            continue
        c = a.obs[col].astype(str).value_counts()
        for lab, n in c.items():
            rows.append({"roi": roi, "label": lab, "n_cells": int(n)})
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def spearman_pairwise(df):
    cols = df.columns
    out = pd.DataFrame(index=cols, columns=cols, dtype=float)
    ranks = df.rank()
    for i, a in enumerate(cols):
        for b in cols[i:]:
            r = ranks[a].corr(ranks[b], method="pearson")
            out.loc[a, b] = r
            out.loc[b, a] = r
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roi-glob", required=True,
                    help="Glob matching ROI output dirs, e.g. 'outputs/autotile_*'")
    ap.add_argument("--out", default="outputs/aggregate", help="output directory")
    args = ap.parse_args()

    roi_dirs = sorted(d for d in glob.glob(args.roi_glob) if Path(d).is_dir())
    if not roi_dirs:
        raise SystemExit(f"no ROI dirs matched {args.roi_glob}")
    print(f"[stage 10] aggregating {len(roi_dirs)} ROIs")
    for d in roi_dirs:
        print(f"  {d}")

    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"; fig_dir.mkdir(exist_ok=True)
    tab_dir = out_dir / "tables"; tab_dir.mkdir(exist_ok=True)

    # 1. Per-ROI IgG neighbour z-scores
    z = collect_igg_rows(roi_dirs)
    if z is None:
        print("[stage 10] no ROIs had stage 6 IgG focal row CSV; nothing to aggregate")
        return
    z.to_csv(tab_dir / "10_cross_roi_IgG_neighbour_zscore.csv")
    print(f"[stage 10] wrote {tab_dir / '10_cross_roi_IgG_neighbour_zscore.csv'}")

    # 2. Per-ROI summaries: mean z, sign consistency, fraction-positive across ROIs
    summary = pd.DataFrame({
        "mean_z":       z.mean(axis=1),
        "median_z":     z.median(axis=1),
        "sign_consistency": np.where(
            z.mean(axis=1) >= 0,
            (z > 0).mean(axis=1),
            (z < 0).mean(axis=1),
        ),
        "n_rois_seen":  (z != 0).sum(axis=1),
    }).sort_values("mean_z", ascending=False)
    summary.to_csv(tab_dir / "10_cross_roi_IgG_neighbour_summary.csv")
    print(f"[stage 10] neighbour summary (top 10 by mean_z):\n"
          f"{summary.head(10).round(3).to_string()}")

    # 3. Pairwise rank correlation of the IgG focal-row z-score pattern
    if z.shape[1] >= 2:
        corr = spearman_pairwise(z)
        corr.to_csv(tab_dir / "10_cross_roi_IgG_pattern_correlation.csv")
        # Triangular heatmap
        fig, ax = plt.subplots(figsize=(8, 7))
        im = ax.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(corr.columns)))
        ax.set_yticks(range(len(corr.index)))
        ax.set_xticklabels(corr.columns, rotation=90, fontsize=7)
        ax.set_yticklabels(corr.index, fontsize=7)
        plt.colorbar(im, ax=ax, label="Spearman rank corr (IgG focal row)")
        ax.set_title("Cross-ROI reproducibility of IgG-neighbour pattern")
        plt.tight_layout()
        plt.savefig(fig_dir / "10_cross_roi_IgG_pattern_correlation.png",
                    dpi=180, bbox_inches="tight")
        plt.close()
        mean_off_diag = (corr.values.sum() - np.trace(corr.values)) / (corr.size - len(corr))
        print(f"[stage 10] mean off-diagonal correlation: {mean_off_diag:.3f}  "
              f"(> 0.5 = strong reproducibility)")

    # 4. Cell-type counts per ROI (long format)
    counts = collect_celltype_counts(roi_dirs)
    if not counts.empty:
        counts.to_csv(tab_dir / "10_cross_roi_celltype_counts.csv", index=False)
        wide = counts.pivot_table(index="label", columns="roi",
                                  values="n_cells", fill_value=0)
        wide["total"] = wide.sum(axis=1)
        wide = wide.sort_values("total", ascending=False)
        wide.to_csv(tab_dir / "10_cross_roi_celltype_counts_wide.csv")
        print(f"[stage 10] cell-type counts (top 10 labels by total):\n"
              f"{wide.head(10).to_string()}")


if __name__ == "__main__":
    main()
