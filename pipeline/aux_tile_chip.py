"""Auxiliary helper: tile a chip into multiple ROI YAML configs.

Generates one configs/<base_name>_NN.yaml per non-empty tile so the chip can
be analysed thoroughly, not as one biased strip.

Reads:  proseg full-chip AnnData (input.proseg_h5ad)
        a "base" config to copy QC/BANKSY/annotation blocks from
Writes: configs/<base_name>_01.yaml, _02.yaml, ...
        configs/<base_name>_tiles_overview.png

Usage:
    python pipeline/aux_tile_chip.py --base-config configs/strip_01_local.yaml \\
        --tile-width 7000 --tile-height 2500 --tile-overlap 500 \\
        --min-cells-per-tile 2000 --base-name autotile

After this writes the YAMLs, run each one:
    for f in configs/autotile_*.yaml; do
        python pipeline/run_all.py --config "$f"
    done

Or in one go later via pipeline/10_aggregate_rois.py.
"""
from pathlib import Path
import argparse
import copy

import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import yaml

import sys; sys.path.insert(0, str(Path(__file__).parent))
from io_utils import load_config


def tile_grid(xy, tile_w, tile_h, overlap, min_cells):
    """Return list of (x_min, x_max, y_min, y_max) tiles covering the chip."""
    x_min_all, x_max_all = float(xy[:, 0].min()), float(xy[:, 0].max())
    y_min_all, y_max_all = float(xy[:, 1].min()), float(xy[:, 1].max())
    step_x = tile_w - overlap
    step_y = tile_h - overlap
    tiles = []
    y = y_min_all
    while y < y_max_all:
        x = x_min_all
        while x < x_max_all:
            x2 = min(x + tile_w, x_max_all)
            y2 = min(y + tile_h, y_max_all)
            mask = ((xy[:, 0] >= x) & (xy[:, 0] < x2) &
                    (xy[:, 1] >= y) & (xy[:, 1] < y2))
            n = int(mask.sum())
            if n >= min_cells:
                tiles.append((int(x), int(x2), int(y), int(y2), n))
            x += step_x
        y += step_y
    return tiles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-config", required=True,
                    help="An existing ROI YAML to copy QC/BANKSY/annotation blocks from.")
    ap.add_argument("--tile-width",  type=int, default=7000, help="tile width (px)")
    ap.add_argument("--tile-height", type=int, default=2500, help="tile height (px)")
    ap.add_argument("--tile-overlap", type=int, default=500, help="overlap between tiles (px)")
    ap.add_argument("--min-cells-per-tile", type=int, default=2000,
                    help="skip tiles with fewer than this many cells")
    ap.add_argument("--base-name", default="autotile",
                    help="prefix for generated config files and roi_id")
    args = ap.parse_args()

    cfg = load_config(args.base_config)
    base_config_path = Path(args.base_config).resolve()
    proseg_path = cfg["input"]["proseg_h5ad"]

    print(f"loading full chip from {proseg_path} ...")
    a = ad.read_h5ad(proseg_path)
    xy = np.asarray(a.obsm["spatial"])
    print(f"chip: {a.n_obs:,} cells, x[{xy[:,0].min():.0f}, {xy[:,0].max():.0f}], "
          f"y[{xy[:,1].min():.0f}, {xy[:,1].max():.0f}]")

    tiles = tile_grid(xy, args.tile_width, args.tile_height,
                      args.tile_overlap, args.min_cells_per_tile)
    print(f"generated {len(tiles)} non-empty tiles")

    configs_dir = base_config_path.parent
    written = []
    for i, (x1, x2, y1, y2, n) in enumerate(tiles, start=1):
        tile_cfg = copy.deepcopy(cfg)
        roi_id = f"{args.base_name}_{i:02d}"
        tile_cfg["roi_id"] = roi_id
        tile_cfg["description"] = (
            f"auto-tiled ROI {i} of {len(tiles)}: "
            f"x[{x1},{x2}] y[{y1},{y2}], ~{n} cells"
        )
        tile_cfg["strip_coords"] = dict(x_min=x1, x_max=x2, y_min=y1, y_max=y2)
        tile_cfg["output_dir"] = f"outputs/{roi_id}"
        # Drop the broad_label_map so stage 4 auto-suggests (clusters differ per tile)
        if "annotation" in tile_cfg and "broad_label_map" in tile_cfg["annotation"]:
            del tile_cfg["annotation"]["broad_label_map"]

        out_path = configs_dir / f"{roi_id}_local.yaml"
        with open(out_path, "w") as f:
            yaml.safe_dump(tile_cfg, f, sort_keys=False)
        written.append(out_path)

    print(f"wrote {len(written)} config files to {configs_dir}/")
    print(f"first: {written[0].name}, last: {written[-1].name}")

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.scatter(xy[:, 0], xy[:, 1], s=0.5, c="grey", alpha=0.25, edgecolors="none")
    for i, (x1, x2, y1, y2, n) in enumerate(tiles, start=1):
        ax.add_patch(Rectangle((x1, y1), x2 - x1, y2 - y1,
                               fill=False, edgecolor="red", linewidth=1.5))
        ax.text(x1, y1 - 150, f"{i:02d} ({n})", color="red", fontsize=8)
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_title(f"chip tiling: {len(tiles)} ROIs covering {a.n_obs:,} cells "
                 f"({args.tile_width}x{args.tile_height} px, overlap {args.tile_overlap})")
    fig_path = configs_dir / f"{args.base_name}_tiles_overview.png"
    plt.tight_layout()
    plt.savefig(fig_path, dpi=140, bbox_inches="tight")
    plt.close()
    print(f"overview figure: {fig_path}")


if __name__ == "__main__":
    main()
