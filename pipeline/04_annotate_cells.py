"""Stage 4: annotate clusters into broad cell types, then fibroblast/vessel/immune subtypes.

Reads: outputs/<roi_id>/adata_banksy.h5ad
       configs/marker_panels.yaml (path from config)
Writes: outputs/<roi_id>/adata_annotated.h5ad
        outputs/<roi_id>/figures/04_*.png
        outputs/<roi_id>/tables/04_*.csv

Substages:
  4a. Broad annotation at canonical_lambda using the broad_label_map in the config.
  4b. Fibroblast subtypes: score the 5 panels on fibroblast/deep_dermal_stromal cells,
      assign highest-scoring panel above the z threshold; otherwise fibroblast_uncertain.
  4c. Vessel labels: score 3 panels on ALL cells, apply if z >= threshold.
  4d. Immune labels: score 4 panels on ALL cells, apply if z >= threshold,
      preserving IgG_producing strict call later (set in stage 5).

Run:
    python pipeline/04_annotate_cells.py --config configs/strip_01.yaml
"""
from pathlib import Path
import argparse
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

import sys; sys.path.insert(0, str(Path(__file__).parent))
from io_utils import load_config, load_marker_panels, out_paths, patch_log1p_base
from stats_utils import score_marker_sets, marker_fraction_table
from plot_utils import spatial_scatter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = out_paths(cfg)
    panels = load_marker_panels(cfg["annotation"]["marker_panels"])
    z_thresh = cfg["annotation"]["min_marker_score_z"]
    raw_map = cfg["annotation"].get("broad_label_map") or {}
    label_map = {str(k): v for k, v in raw_map.items()}
    banksy_key = f"banksy_l{cfg['canonical_lambda']}"

    print(f"[stage 4] roi_id      : {cfg['roi_id']}")
    print(f"[stage 4] using BANKSY: {banksy_key}")
    print(f"[stage 4] z threshold : {z_thresh}")

    a = ad.read_h5ad(paths["banksy"])
    patch_log1p_base(a)   # scanpy 1.9 h5ad round-trip workaround
    print(f"[stage 4] loaded      : {a.shape}")

    # ---- 4a. Broad annotation ---------------------------------------------
    if banksy_key not in a.obs.columns:
        raise KeyError(f"{banksy_key} not in adata.obs. Re-run stage 3 with that lambda.")

    # Marker table per BANKSY cluster (used both for sanity-check and for
    # auto-suggesting labels when broad_label_map is empty).
    sc.tl.rank_genes_groups(a, groupby=banksy_key, method="wilcoxon", use_raw=False)
    top_genes = pd.DataFrame({
        c: list(a.uns["rank_genes_groups"]["names"][c][:25])
        for c in a.uns["rank_genes_groups"]["names"].dtype.names
    })
    top_genes.to_csv(paths["tab"] / "04a_top_markers_per_cluster.csv", index=False)
    print(f"[4a] wrote top markers per cluster")

    # If the YAML provided a label map, use it. Otherwise auto-suggest from
    # marker overlap with the broad_markers panels and write suggestions for
    # the user to review.
    cluster_ids = list(a.obs[banksy_key].astype(str).unique())
    missing = [c for c in cluster_ids if c not in label_map]
    if missing:
        print(f"[4a] {len(missing)} cluster(s) not in broad_label_map; auto-suggesting from markers")
        broad_panels = panels.get("broad_markers", {})
        suggestions = {}
        rows = []
        for cl in cluster_ids:
            cl_top = set(top_genes[cl].head(25))
            best_panel, best_overlap, best_n = None, 0, 0
            for panel_name, panel_genes in broad_panels.items():
                overlap = cl_top.intersection(panel_genes)
                if len(overlap) > best_overlap:
                    best_panel, best_overlap, best_n = panel_name, len(overlap), len(panel_genes)
            suggestions[cl] = best_panel or "unassigned"
            rows.append({
                "cluster": cl,
                "n_cells": int((a.obs[banksy_key].astype(str) == cl).sum()),
                "suggested_label": suggestions[cl],
                "marker_overlap": best_overlap,
                "panel_size": best_n,
                "top_markers": ",".join(top_genes[cl].head(10)),
            })
        sug_df = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
        sug_df.to_csv(paths["tab"] / "04a_broad_label_suggestions.csv", index=False)
        print(f"[4a] auto-suggested labels:\n{sug_df.to_string(index=False)}")
        print(f"[4a] To override, add a broad_label_map block to your YAML and re-run stage 4.")
        for cl, lab in suggestions.items():
            label_map.setdefault(cl, lab)

    a.obs["celltype"] = a.obs[banksy_key].astype(str).map(label_map)
    n_unmapped = int(a.obs["celltype"].isna().sum())
    if n_unmapped > 0:
        raise ValueError(
            f"{n_unmapped} cells still unmapped after auto-suggestion. "
            f"Cluster -> count summary:\n{a.obs[banksy_key].value_counts()}"
        )
    print(f"[4a] broad cell counts:\n{a.obs['celltype'].value_counts().to_string()}")

    spatial_scatter(
        a, "celltype",
        title=f"{cfg['roi_id']}: broad cell type annotation",
        filename=paths["fig"] / "04a_broad_celltype_spatial.png",
    )

    # ---- 4b. Fibroblast subtypes -----------------------------------------
    fib_mask = a.obs["celltype"].isin(["fibroblast", "deep_dermal_stromal"]).values
    fib = a[fib_mask].copy()
    print(f"[4b] scoring {fib.n_obs:,} fibroblast/stromal cells")
    fib_panels = panels["fibroblast_subtype_markers"]
    fib_info = score_marker_sets(fib, fib_panels, prefix="fibro")
    fib_info.to_csv(paths["tab"] / "04b_fibroblast_marker_presence.csv", index=False)
    print(f"[4b] fibroblast panel coverage:\n{fib_info[['marker_set', 'n_present', 'n_missing']].to_string(index=False)}")

    score_cols = [f"fibro_{k}" for k in fib_panels]
    scores = fib.obs[score_cols].copy()
    scores.columns = [c.replace("fibro_", "") for c in scores.columns]
    fib.obs["fib_subtype"] = scores.idxmax(axis=1)
    fib.obs.loc[scores.max(axis=1) < z_thresh, "fib_subtype"] = "fibroblast_uncertain"
    print(f"[4b] fibroblast subtype counts:\n{fib.obs['fib_subtype'].value_counts().to_string()}")

    # Push subtype labels back onto the main AnnData
    a.obs["celltype_detailed_v0"] = a.obs["celltype"].astype(str)
    a.obs.loc[fib.obs.index, "celltype_detailed_v0"] = fib.obs["fib_subtype"].values

    spatial_scatter(
        a, "celltype_detailed_v0",
        title=f"{cfg['roi_id']}: fibroblast subtype proposal",
        filename=paths["fig"] / "04b_fibroblast_subtype_spatial.png",
    )

    # ---- 4c. Vessel labels -----------------------------------------------
    vasc_info = score_marker_sets(a, panels["vessel_markers"], prefix="vascular")
    vasc_info.to_csv(paths["tab"] / "04c_vessel_marker_presence.csv", index=False)
    print(f"[4c] vessel panel coverage:\n{vasc_info[['marker_set', 'n_present', 'n_missing']].to_string(index=False)}")

    # ---- 4d. Immune labels + final detailed label vector -----------------
    imm_info = score_marker_sets(a, panels["immune_markers"], prefix="immune")
    imm_info.to_csv(paths["tab"] / "04d_immune_marker_presence.csv", index=False)
    print(f"[4d] immune panel coverage:\n{imm_info[['marker_set', 'n_present', 'n_missing']].to_string(index=False)}")

    candidate = a.obs["celltype_detailed_v0"].astype(str).copy()
    # Priority order: highest-specificity panels first
    override_rules = [
        ("immune_mast_cell",                 "mast_cell_candidate"),
        ("immune_t_cell",                    "t_cell_candidate"),
        ("immune_macrophage",                "macrophage_candidate"),
        ("immune_b_cell",                    "b_cell_candidate"),
        ("vascular_endothelial",             "blood_endothelial_candidate"),
        ("vascular_lymphatic_endothelial",   "lymphatic_endothelial_candidate"),
        ("vascular_pericyte_smooth_muscle",  "pericyte_smooth_muscle_candidate"),
    ]
    for sc_col, lab in override_rules:
        if sc_col in a.obs.columns and a.obs[sc_col].notna().any():
            m = a.obs[sc_col].values >= z_thresh
            candidate.loc[m] = lab
            print(f"[4d]   {lab:<32} <- {int(m.sum())} cells")

    a.obs["celltype_detailed_v1"] = pd.Categorical(candidate)
    print(f"[4d] celltype_detailed_v1 counts:\n{a.obs['celltype_detailed_v1'].value_counts().to_string()}")

    # Marker fraction sanity check
    selected = [
        "POSTN", "ASPN", "CXCL14", "PI16", "DPT", "ACTA2", "TAGLN",
        "PECAM1", "VWF", "LYVE1", "PROX1", "RGS5", "PDGFRB",
        "LYZ", "CD68", "CD163", "TPSAB1", "CD3D", "MS4A1", "CD79A",
        "IGHG1", "IGHG4", "IGKC", "MZB1",
    ]
    frac = marker_fraction_table(a, "celltype_detailed_v1", selected)
    frac.to_csv(paths["tab"] / "04d_detailed_label_marker_fraction.csv")
    print(f"[4d] wrote per-label marker fraction table")

    spatial_scatter(
        a, "celltype_detailed_v1",
        title=f"{cfg['roi_id']}: detailed cell-type annotation (pre-IgG)",
        filename=paths["fig"] / "04d_detailed_celltype_v1_spatial.png",
    )

    a.write_h5ad(paths["annotated"])
    print(f"[stage 4] wrote: {paths['annotated']}")


if __name__ == "__main__":
    main()
