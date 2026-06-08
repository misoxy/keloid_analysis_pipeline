"""Stage 5: identify IgG-rich candidate cells with multi-definition sensitivity.

Reads:  outputs/<roi_id>/adata_annotated.h5ad
Writes: outputs/<roi_id>/adata_with_igg.h5ad
        outputs/<roi_id>/tables/05_igg_*.csv
        outputs/<roi_id>/figures/05_igg_*.png

Logic:
  1. Sub-cluster the deep_dermal_stromal compartment at higher resolution.
  2. Pick the sub-cluster with the highest IGHG1+ fraction as the IgG-rich
     candidate sub-cluster (Set B).
  3. Define alternative IgG sets (A: IGHG1>0; C: top 2.5% Ig score;
     D = B AND IGHG1+ = strict).
  4. Marker fractions per set so the lab member sees whether heavy / light
     chain enrichment is real.
  5. Apply the honest "IgG_rich_candidate" label = Set D into celltype_hc.
  6. Apply Set D into celltype_detailed_v1 (preserving the IgG cluster
     against vessel/immune overrides from stage 4).

Run:
    python pipeline/05_igg_detection.py --config configs/strip_01.yaml
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
from io_utils import load_config, out_paths, patch_log1p_base
from stats_utils import gene_vector
from plot_utils import spatial_scatter


PLASMA_LABEL = "IgG_rich_candidate"   # honest, not "plasma_cell"


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
    ig_cfg = cfg["igg"]
    SEED = 0

    print(f"[stage 5] roi_id: {cfg['roi_id']}")
    a = ad.read_h5ad(paths["annotated"])
    patch_log1p_base(a)   # scanpy 1.9 h5ad round-trip workaround
    print(f"[stage 5] loaded: {a.shape}")

    # ---- 1. Sub-cluster deep_dermal_stromal ------------------------------
    sub = a[a.obs["celltype"] == "deep_dermal_stromal"].copy()
    print(f"[5.1] sub-clustering {sub.n_obs:,} deep_dermal_stromal cells")
    sc.pp.highly_variable_genes(sub, n_top_genes=ig_cfg["sub_n_hvg"], flavor="seurat", subset=False)
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=ig_cfg["sub_n_pcs"])
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=ig_cfg["sub_n_pcs"])
    sc.tl.leiden(sub, resolution=ig_cfg["sub_leiden_resolution"],
                 random_state=SEED, key_added="sub_leiden")

    ighg1 = gene_vector(sub, "IGHG1")
    sub.obs["ighg1_pos"] = ighg1 > 0
    enrich = sub.obs.groupby("sub_leiden")["ighg1_pos"].agg(["sum", "count"])
    enrich["fraction_ighg1_pos"] = enrich["sum"] / enrich["count"]
    enrich.to_csv(paths["tab"] / "05_subcluster_ighg1_enrichment.csv")
    print(f"[5.1] IGHG1+ fraction per sub-cluster:\n{enrich.sort_values('fraction_ighg1_pos', ascending=False).to_string()}")

    best_sub = enrich["fraction_ighg1_pos"].idxmax()
    best_frac = float(enrich.loc[best_sub, "fraction_ighg1_pos"])
    best_n = int(enrich.loc[best_sub, "count"])
    print(f"[5.1] selected sub-cluster {best_sub}: {best_n} cells, {100*best_frac:.1f}% IGHG1+")
    igg_idx = sub.obs.index[sub.obs["sub_leiden"] == best_sub]

    # ---- 2. Define candidate sets ----------------------------------------
    ighg1_main = gene_vector(a, "IGHG1")
    set_A = ighg1_main > 0
    set_B = a.obs.index.isin(igg_idx)

    panel_in = [g for g in ig_cfg["ig_score_panel"] if g in a.var_names]
    sc.tl.score_genes(a, gene_list=panel_in, score_name="ig_score", random_state=SEED)
    set_C = (a.obs["ig_score"].values > a.obs["ig_score"].quantile(0.975))
    set_D = set_A & set_B

    sets = {
        "A (IGHG1>0)":          set_A,
        "B (sub-cluster)":      set_B,
        "C (Ig score top 2.5%)": set_C,
        "D (B AND IGHG1>0)":    set_D,
        "background":           ~(set_A | set_B | set_C | set_D),
    }
    sizes = pd.Series({k: int(v.sum()) for k, v in sets.items()}, name="n")
    print(f"[5.2] candidate set sizes:\n{sizes.to_string()}")
    sizes.to_csv(paths["tab"] / "05_igg_set_sizes.csv", header=True)

    # ---- 3. Marker validation per set ------------------------------------
    panel = ["IGHG1","IGHG2","IGHG3","IGHG4","IGKC","IGLC1","IGLC2","IGLC3",
             "JCHAIN","MZB1","XBP1","SDC1","CD79A","CD38"]
    frac = pd.DataFrame({
        name: fraction_expressing(a, mask, panel) for name, mask in sets.items()
    })
    frac.to_csv(paths["tab"] / "05_igg_set_marker_fraction.csv")
    print(f"[5.3] marker fraction per set:\n{frac.round(3).to_string()}")

    # ---- 4. Apply honest label -------------------------------------------
    HC_SET = set_D if set_D.sum() >= 30 else set_B
    print(f"[5.4] using {'strict Set D' if HC_SET is set_D else 'fallback Set B'}: {int(HC_SET.sum())} cells")

    a.obs["celltype_hc"] = a.obs["celltype"].astype(str)
    a.obs.loc[HC_SET, "celltype_hc"] = PLASMA_LABEL
    a.obs["ighg1_pos"] = set_A   # keep loose evidence available

    # Inject into celltype_detailed_v1 so downstream stages see the IgG label
    cdv1 = a.obs["celltype_detailed_v1"].astype(str).copy()
    cdv1.loc[HC_SET] = PLASMA_LABEL
    a.obs["celltype_detailed_v1"] = pd.Categorical(cdv1)
    print(f"[5.4] final celltype_hc:\n{a.obs['celltype_hc'].value_counts().to_string()}")

    # ---- 5. IgG sensitivity: alternative definitions --------------------
    def make_set(name):
        if name == "strict":
            return set_D
        if name == "ig_score_top_2p5pct":
            return set_C
        if name == "ighg1_and_ighg4":
            v4 = gene_vector(a, "IGHG4")
            return set_A & (v4 > 0)
        if name == "ighg1_or_ighg4_plus_igkc":
            v4 = gene_vector(a, "IGHG4"); vk = gene_vector(a, "IGKC")
            return ((ighg1_main > 0) | (v4 > 0)) & (vk > 0)
        raise ValueError(f"unknown IgG definition: {name}")

    sens_rows = []
    for name in ig_cfg["sensitivity_definitions"]:
        s = make_set(name)
        sens_rows.append({"definition": name, "n_cells": int(s.sum())})
    sens_df = pd.DataFrame(sens_rows)
    sens_df.to_csv(paths["tab"] / "05_igg_sensitivity_definitions.csv", index=False)
    print(f"[5.5] sensitivity definitions:\n{sens_df.to_string(index=False)}")

    # Spatial highlight figure
    spatial_scatter(
        a, "celltype_hc",
        title=f"{cfg['roi_id']}: broad annotation with {PLASMA_LABEL} highlighted",
        filename=paths["fig"] / "05_igg_spatial.png",
        highlight=[PLASMA_LABEL],
    )

    a.write_h5ad(paths["with_igg"])
    print(f"[stage 5] wrote: {paths['with_igg']}")


if __name__ == "__main__":
    main()
