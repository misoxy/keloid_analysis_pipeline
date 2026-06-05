"""QC, gene filter, normalise, and pick highly variable genes for one strip.

Reads outputs/<roi_id>/adata_strip.h5ad, applies the QC settings from the YAML
config, and writes outputs/<roi_id>/adata_normalised.h5ad.

We use the same conventions as the existing proseg pipeline:
  - Cap UMI per cell at 99th percentile (drops overlap doublets).
  - Drop near-empty cells and barely-seen genes.
  - Pick 3000 HVGs with seurat_v3, then drop structural gene families
    (COL, KRT, MT-, RP, HB) so they don't hijack downstream clustering.
  - Normalise to a target total + log1p. We do NOT scale here -- BANKSY
    handles its own scaling internally.

Run:
    python pipeline/02_qc_normalise.py --config configs/strip_01.yaml
"""

from pathlib import Path
import argparse
import yaml

import numpy as np
import scanpy as sc
import anndata as ad


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    out_dir = Path(cfg["output_dir"])
    in_path  = out_dir / "adata_strip.h5ad"
    out_path = out_dir / "adata_normalised.h5ad"

    qc = cfg["qc"]

    print(f"roi_id : {cfg['roi_id']}")
    print(f"input  : {in_path}")

    a = ad.read_h5ad(in_path)
    a.var_names_make_unique()
    print(f"loaded : {a.n_obs:,} cells x {a.n_vars:,} genes")

    # Keep a raw counts layer so HVG can compute on raw counts later.
    a.layers["counts"] = a.X.copy()

    # 1. Drop the top 1% UMI cells (likely doublets / overlapping cells)
    umi = np.asarray(a.X.sum(1)).flatten()
    cap = np.percentile(umi, qc["cap_umi_percentile"])
    keep = umi <= cap
    a = a[keep].copy()
    print(f"  cap UMI at p{qc['cap_umi_percentile']} = {int(cap):,} -> {a.n_obs:,} cells")

    # 2. Drop near-empty cells and barely-seen genes
    sc.pp.filter_cells(a, min_counts=qc["min_counts_per_cell"])
    sc.pp.filter_genes(a, min_cells=qc["min_cells_per_gene"])
    print(f"  after cell/gene filter: {a.n_obs:,} cells x {a.n_vars:,} genes")

    # 3. Pick highly variable genes on raw counts using seurat_v3
    sc.pp.highly_variable_genes(
        a, n_top_genes=qc["n_top_hvg"],
        flavor="seurat_v3", layer="counts",
    )

    # 4. Force structural families out of the HVG set so they don't dominate
    struct = a.var_names.str.contains(qc["structural_regex"], regex=True)
    n_dropped = int((a.var["highly_variable"] & struct).sum())
    a.var.loc[struct, "highly_variable"] = False
    n_hvg = int(a.var["highly_variable"].sum())
    print(f"  HVG: {n_hvg} kept, {n_dropped} structural dropped")

    # 5. Normalise to a common total + log1p. NO scaling here.
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    # 6. Add explicit x and y columns BANKSY's initialize_banksy expects.
    a.obs["x"] = a.obsm["spatial"][:, 0]
    a.obs["y"] = a.obsm["spatial"][:, 1]

    out_dir.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(out_path)
    print(f"wrote  : {out_path}  ({out_path.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
