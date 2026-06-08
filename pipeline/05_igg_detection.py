"""Stage 5: generic focal-aggregate detection (default: IgG in keloid).

Reads:  outputs/<roi_id>/adata_annotated.h5ad
Writes: outputs/<roi_id>/adata_with_igg.h5ad
        outputs/<roi_id>/tables/05_*.csv
        outputs/<roi_id>/figures/05_*.png

Generic logic, controlled entirely by the `focal:` block in the YAML config:

  1. Sub-cluster the target_compartment (a broad celltype) at higher resolution.
  2. Pick the sub-cluster with the highest fraction of primary_marker-positive
     cells. This is "Set B" (sub-cluster).
  3. Define candidate sets relative to the primary_marker:
       A: primary_marker > 0 anywhere on the strip
       B: sub-cluster from step 2
       C: top 2.5% on the composite ig_score_panel
       D: A AND B  -> the strict definition
  4. Compute marker fractions per set using validation_panel.
  5. Apply target_label to the high-confidence set (D if >= 30 cells, else B).
  6. Sensitivity analysis across alternative definitions named in
     sensitivity_definitions.

Defaults to IgG detection in deep_dermal_stromal for the keloid project, but
every keloid-specific piece is configurable. Set focal.enabled: false to skip.

Run:
    python pipeline/05_igg_detection.py --config configs/strip_01.yaml
"""
from pathlib import Path
import argparse
import shutil
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

import sys; sys.path.insert(0, str(Path(__file__).parent))
from io_utils import load_config, out_paths, patch_log1p_base
from stats_utils import gene_vector
from plot_utils import spatial_scatter


def fraction_expressing(adata, mask, genes):
    out = {}
    for g in genes:
        v = gene_vector(adata, g)
        out[g] = float((v[mask] > 0).mean()) if (v is not None and mask.sum() > 0) else np.nan
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = out_paths(cfg)
    # Back-compat: read either "focal" or legacy "igg" block.
    fcfg = cfg.get("focal") or cfg.get("igg") or {}
    SEED = 0

    enabled = fcfg.get("enabled", True)
    if not enabled:
        print(f"[stage 5] focal.enabled = false; skipping. Copying annotated -> with_igg.")
        shutil.copy(paths["annotated"], paths["with_igg"])
        return

    target_compartment = fcfg.get("target_compartment", "deep_dermal_stromal")
    target_label       = fcfg.get("target_label", "IgG_rich_candidate")
    primary_marker     = fcfg.get("primary_marker", "IGHG1")
    secondary_marker   = fcfg.get("secondary_marker", "IGHG4")
    light_marker       = fcfg.get("light_chain_marker", "IGKC")
    validation_panel   = fcfg.get("validation_panel",
                                  ["IGHG1","IGHG2","IGHG3","IGHG4","IGKC","IGLC1","IGLC2","IGLC3",
                                   "JCHAIN","MZB1","XBP1","SDC1","CD79A","CD38"])
    score_panel        = fcfg.get("ig_score_panel",
                                  ["IGHG1","IGHG2","IGHG3","IGHG4","IGKC","JCHAIN"])
    sens_defs          = fcfg.get("sensitivity_definitions",
                                  ["strict", "score_top_2p5pct",
                                   "primary_and_secondary", "primary_or_secondary_plus_light"])

    print(f"[stage 5] roi_id           : {cfg['roi_id']}")
    print(f"[stage 5] target_compartment: {target_compartment}")
    print(f"[stage 5] target_label     : {target_label}")
    print(f"[stage 5] primary_marker   : {primary_marker}")

    a = ad.read_h5ad(paths["annotated"])
    patch_log1p_base(a)
    print(f"[stage 5] loaded          : {a.shape}")

    # ---- 1. Sub-cluster target_compartment -------------------------------
    if target_compartment not in a.obs["celltype"].astype(str).unique():
        print(f"[stage 5] WARNING: '{target_compartment}' not in broad celltype labels. Skipping.")
        a.write_h5ad(paths["with_igg"])
        return

    sub = a[a.obs["celltype"] == target_compartment].copy()
    print(f"[5.1] sub-clustering {sub.n_obs:,} '{target_compartment}' cells")
    sc.pp.highly_variable_genes(sub, n_top_genes=fcfg.get("sub_n_hvg", 2000),
                                flavor="seurat", subset=False)
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=fcfg.get("sub_n_pcs", 30))
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=fcfg.get("sub_n_pcs", 30))
    sc.tl.leiden(sub, resolution=fcfg.get("sub_leiden_resolution", 0.6),
                 random_state=SEED, key_added="sub_leiden")

    primary_vec_sub = gene_vector(sub, primary_marker)
    if primary_vec_sub is None:
        print(f"[stage 5] WARNING: primary_marker '{primary_marker}' not in var_names. Skipping.")
        a.write_h5ad(paths["with_igg"])
        return
    sub.obs["primary_pos"] = primary_vec_sub > 0
    enrich = sub.obs.groupby("sub_leiden")["primary_pos"].agg(["sum", "count"])
    enrich["fraction_primary_pos"] = enrich["sum"] / enrich["count"]
    enrich.to_csv(paths["tab"] / "05_subcluster_primary_enrichment.csv")
    print(f"[5.1] {primary_marker}+ fraction per sub-cluster:\n"
          f"{enrich.sort_values('fraction_primary_pos', ascending=False).to_string()}")

    best_sub = enrich["fraction_primary_pos"].idxmax()
    best_frac = float(enrich.loc[best_sub, "fraction_primary_pos"])
    best_n = int(enrich.loc[best_sub, "count"])
    print(f"[5.1] selected sub-cluster {best_sub}: {best_n} cells, "
          f"{100*best_frac:.1f}% {primary_marker}+")
    target_idx = sub.obs.index[sub.obs["sub_leiden"] == best_sub]

    # ---- 2. Build candidate sets -----------------------------------------
    primary_vec = gene_vector(a, primary_marker)
    set_A = primary_vec > 0
    set_B = a.obs.index.isin(target_idx)

    panel_in = [g for g in score_panel if g in a.var_names]
    sc.tl.score_genes(a, gene_list=panel_in, score_name="composite_score", random_state=SEED)
    set_C = (a.obs["composite_score"].values > a.obs["composite_score"].quantile(0.975))
    set_D = set_A & set_B

    sets = {
        f"A ({primary_marker}>0)":             set_A,
        "B (sub-cluster)":                     set_B,
        "C (composite score top 2.5%)":        set_C,
        f"D (B AND {primary_marker}>0)":       set_D,
        "background":                          ~(set_A | set_B | set_C | set_D),
    }
    sizes = pd.Series({k: int(v.sum()) for k, v in sets.items()}, name="n")
    print(f"[5.2] candidate set sizes:\n{sizes.to_string()}")
    sizes.to_csv(paths["tab"] / "05_set_sizes.csv", header=True)

    # ---- 3. Marker validation per set ------------------------------------
    frac = pd.DataFrame({
        name: fraction_expressing(a, mask, validation_panel)
        for name, mask in sets.items()
    })
    frac.to_csv(paths["tab"] / "05_set_marker_fraction.csv")
    print(f"[5.3] marker fraction per set:\n{frac.round(3).to_string()}")

    # ---- 4. Apply honest label -------------------------------------------
    HC_SET = set_D if set_D.sum() >= 30 else set_B
    set_name = "strict Set D" if HC_SET is set_D else "fallback Set B"
    print(f"[5.4] applying '{target_label}' to {set_name}: {int(HC_SET.sum())} cells")

    a.obs["celltype_hc"] = a.obs["celltype"].astype(str)
    a.obs.loc[HC_SET, "celltype_hc"] = target_label
    a.obs["primary_pos"] = set_A

    cdv1 = a.obs["celltype_detailed_v1"].astype(str).copy()
    cdv1.loc[HC_SET] = target_label
    a.obs["celltype_detailed_v1"] = pd.Categorical(cdv1)
    print(f"[5.4] final celltype_hc:\n{a.obs['celltype_hc'].value_counts().to_string()}")

    # ---- 5. Sensitivity definitions --------------------------------------
    def make_set(name):
        if name == "strict":
            return set_D
        if name == "score_top_2p5pct":
            return set_C
        if name == "primary_and_secondary":
            v2 = gene_vector(a, secondary_marker)
            if v2 is None: return set_A.copy()
            return set_A & (v2 > 0)
        if name == "primary_or_secondary_plus_light":
            v2 = gene_vector(a, secondary_marker); vL = gene_vector(a, light_marker)
            if v2 is None or vL is None: return set_A.copy()
            return ((primary_vec > 0) | (v2 > 0)) & (vL > 0)
        # Legacy aliases (back-compat with strip_01.yaml prior wording)
        if name == "ig_score_top_2p5pct":
            return set_C
        if name == "ighg1_and_ighg4":
            v2 = gene_vector(a, "IGHG4")
            if v2 is None: return set_A.copy()
            return set_A & (v2 > 0)
        if name == "ighg1_or_ighg4_plus_igkc":
            v2 = gene_vector(a, "IGHG4"); vL = gene_vector(a, "IGKC")
            if v2 is None or vL is None: return set_A.copy()
            return ((primary_vec > 0) | (v2 > 0)) & (vL > 0)
        raise ValueError(f"unknown sensitivity definition: {name}")

    sens_rows = []
    for name in sens_defs:
        s = make_set(name)
        sens_rows.append({"definition": name, "n_cells": int(s.sum())})
    sens_df = pd.DataFrame(sens_rows)
    sens_df.to_csv(paths["tab"] / "05_sensitivity_definitions.csv", index=False)
    print(f"[5.5] sensitivity definitions:\n{sens_df.to_string(index=False)}")

    spatial_scatter(
        a, "celltype_hc",
        title=f"{cfg['roi_id']}: broad annotation with {target_label} highlighted",
        filename=paths["fig"] / "05_focal_spatial.png",
        highlight=[target_label],
    )

    a.write_h5ad(paths["with_igg"])
    print(f"[stage 5] wrote: {paths['with_igg']}")


if __name__ == "__main__":
    main()
