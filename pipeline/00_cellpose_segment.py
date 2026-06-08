"""Stage 0 (upstream, server-intended): Cellpose nuclei segmentation on H&E.

Reads:  H&E TIFF path from config (upstream.he_tif)
Writes: outputs/<chip_id>/nuclei_mask.npy.gz  (uint32 label image, full chip resolution)
        outputs/<chip_id>/figures/00_cellpose_qc.png

WARNING: this is a SLOW stage. On a full chip H&E (~20k x 23k pixels) it takes
hours on a laptop CPU. Designed to run on the server with GPU or in tiled
mode. For laptop dev, use --tile-size 1024 to run on a small patch only.

Cellpose v3 nuclei model is the default. The cpsam model is ~100x slower on
CPU and not worth it for iteration.

Run:
    # Full chip (server, intended)
    python pipeline/00_cellpose_segment.py --config configs/strip_01.yaml

    # Laptop dev: one 1024x1024 patch at the centre
    python pipeline/00_cellpose_segment.py --config configs/strip_01.yaml --tile-size 1024
"""
from pathlib import Path
import argparse
import gzip
import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import tifffile
import matplotlib.pyplot as plt

import sys; sys.path.insert(0, str(Path(__file__).parent))
from io_utils import load_config


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--tile-size", type=int, default=None,
                        help="Optional: run on a centre tile of this size (laptop dev).")
    parser.add_argument("--diameter", type=float, default=None,
                        help="Expected nuclear diameter in pixels. Auto-estimated if not given.")
    args = parser.parse_args()

    cfg = load_config(args.config)
    up = cfg.get("upstream", {})
    he_tif = Path(up.get("he_tif", ""))
    chip_id = up.get("chip_id", cfg["roi_id"])
    out_dir = Path(up.get("upstream_out_dir", f"outputs/{chip_id}_upstream"))
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"; fig_dir.mkdir(exist_ok=True)
    out_mask = out_dir / "nuclei_mask.npy.gz"

    print(f"[stage 0] chip_id   : {chip_id}")
    print(f"[stage 0] H&E input : {he_tif}")
    if not he_tif.exists():
        sys.exit(f"H&E TIFF not found at {he_tif}. Set upstream.he_tif in the config.")

    print(f"[stage 0] loading H&E (this can take a minute)...")
    t0 = time.time()
    he = tifffile.imread(str(he_tif))
    print(f"[stage 0] H&E shape: {he.shape} dtype={he.dtype} ({(time.time()-t0):.1f}s)")

    # Optional centre tile for laptop dev
    if args.tile_size:
        cx, cy = he.shape[1] // 2, he.shape[0] // 2
        h = args.tile_size // 2
        he = he[cy-h:cy+h, cx-h:cx+h]
        print(f"[stage 0] tile mode: {he.shape}")

    # If RGB, take greyscale of the green channel (nuclei look best in green inverted)
    if he.ndim == 3:
        img = he[..., 1].astype(np.float32)
    else:
        img = he.astype(np.float32)

    print(f"[stage 0] loading Cellpose...")
    from cellpose import models
    # Cellpose v3 nuclei model. Falls back to CPU if no GPU (laptop case).
    model = models.Cellpose(gpu=False, model_type="nuclei")

    diam = args.diameter or up.get("cellpose_diameter", None)
    print(f"[stage 0] running Cellpose nuclei segmentation "
          f"(diameter={diam if diam else 'auto'})...")
    t1 = time.time()
    masks, flows, styles, diams = model.eval(
        img, diameter=diam, channels=[0, 0],
        flow_threshold=up.get("cellpose_flow_threshold", 0.4),
        cellprob_threshold=up.get("cellpose_cellprob_threshold", 0.0),
    )
    print(f"[stage 0] Cellpose done in {(time.time()-t1):.1f}s")
    print(f"[stage 0] estimated diameter : {diams}")
    print(f"[stage 0] nuclei detected    : {int(masks.max())}")

    # Save mask
    print(f"[stage 0] writing {out_mask} ...")
    with gzip.open(out_mask, "wb", compresslevel=1) as f:
        np.save(f, masks.astype(np.uint32))
    print(f"[stage 0] mask file size: {out_mask.stat().st_size/1e6:.1f} MB")

    # QC figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    axes[0].imshow(img, cmap="gray"); axes[0].set_title("H&E (green channel)")
    axes[0].axis("off")
    axes[1].imshow(masks, cmap="prism"); axes[1].set_title(f"Cellpose mask: {int(masks.max())} nuclei")
    axes[1].axis("off")
    plt.tight_layout()
    plt.savefig(fig_dir / "00_cellpose_qc.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[stage 0] QC figure: {fig_dir / '00_cellpose_qc.png'}")

    print(f"[stage 0] total runtime: {(time.time()-t0):.1f}s")


if __name__ == "__main__":
    main()
