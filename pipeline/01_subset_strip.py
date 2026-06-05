"""Subset the full-chip proseg AnnData to one tissue strip.

Reads a YAML config that says where the strip is (x_min, x_max, y_min, y_max
in chip pixel coordinates) and which input file to use, then writes the
subset to outputs/<roi_id>/adata_strip.h5ad.

Run:
    python pipeline/01_subset_strip.py --config configs/strip_01.yaml
"""

from pathlib import Path
import argparse
import yaml

import anndata as ad
import numpy as np


def load_config(path):
    # Read the YAML config and return it as a plain dict.
    with open(path) as f:
        return yaml.safe_load(f)


def subset_to_strip(adata, x_min, x_max, y_min, y_max):
    # AnnData.obsm['spatial'] is an (n_cells, 2) array of (x, y) chip coords.
    # We keep only cells whose centre lies inside the strip rectangle.
    xy = adata.obsm["spatial"]
    in_strip = (
        (xy[:, 0] >= x_min) & (xy[:, 0] <= x_max) &
        (xy[:, 1] >= y_min) & (xy[:, 1] <= y_max)
    )
    n_in = int(in_strip.sum())
    n_out = int((~in_strip).sum())
    print(f"  cells inside strip : {n_in:,}")
    print(f"  cells outside strip: {n_out:,}")
    return adata[in_strip].copy()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True,
                        help="path to a strip YAML config")
    args = parser.parse_args()

    cfg = load_config(args.config)
    print(f"roi_id  : {cfg['roi_id']}")
    print(f"input   : {cfg['input']['proseg_h5ad']}")

    sc = cfg["strip_coords"]
    print(f"strip   : x=[{sc['x_min']}, {sc['x_max']}]  "
          f"y=[{sc['y_min']}, {sc['y_max']}]")

    # Load the full-chip proseg AnnData.
    a = ad.read_h5ad(cfg["input"]["proseg_h5ad"])
    print(f"loaded  : {a.n_obs:,} cells x {a.n_vars:,} genes")

    # Subset to the strip.
    a_strip = subset_to_strip(
        a,
        x_min=sc["x_min"], x_max=sc["x_max"],
        y_min=sc["y_min"], y_max=sc["y_max"],
    )

    # Write the subset.
    out_dir = Path(cfg["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "adata_strip.h5ad"
    a_strip.write_h5ad(out_path)
    print(f"wrote   : {out_path}  ({out_path.stat().st_size/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
