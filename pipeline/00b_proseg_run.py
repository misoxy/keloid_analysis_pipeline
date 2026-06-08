"""Stage 0b (upstream, server-intended): build transcripts CSV + run proseg + load to AnnData.

Reads:  outputs/<chip_id>_upstream/nuclei_mask.npy.gz  (from stage 0)
        Stereo-seq tissue .gef file (from config upstream.tissue_gef)
        proseg binary on PATH or at upstream.proseg_binary
Writes: outputs/<chip_id>_upstream/transcripts.csv.gz
        outputs/<chip_id>_upstream/output/proseg-*.{csv,mtx}.gz
        outputs/<chip_id>_upstream/proseg_full_raw.h5ad   <- the file all downstream stages read

WARNING: this is the slowest pipeline step.
  - building transcripts.csv: 5-10 min on laptop (streams gene-by-gene from .gef)
  - running proseg: 30-90 min depending on chip
  - loading to AnnData: 2-5 min
Intended for server runs. For laptop dev, use a small chip subset.

Critical proseg gotchas (encoded in this script so you don't re-discover them):
  - --output-path is PREPENDED to every output filename. Filenames must be bare.
  - cell_id MUST come from the mask (0 = unassigned). Don't use raw .gef IDs.
  - overlaps_nucleus MUST come from the mask. Don't use the .gef column.
  - --burnin-samples 100 --samples 100 is the validated convergence setting.

Run:
    python pipeline/00b_proseg_run.py --config configs/strip_01.yaml
"""
from pathlib import Path
import argparse
import gzip
import subprocess
import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import h5py

import sys; sys.path.insert(0, str(Path(__file__).parent))
from io_utils import load_config


def build_transcripts_csv(tissue_gef, mask, out_csv, chunk_genes=200):
    """Stream tissue.gef gene-by-gene, derive cell_id + overlaps_nucleus from the mask,
    write the Xenium-style CSV proseg expects.
    """
    H, W = mask.shape
    print(f"[stage 0b] mask shape: ({H}, {W})  with {int(mask.max())} cells")

    with h5py.File(tissue_gef, "r") as f:
        exp = f["geneExp/bin1/expression"]
        gene = f["geneExp/bin1/gene"]
        total = exp.shape[0]
        print(f"[stage 0b] tissue.gef has {total:,} transcript records")
        gene_names = [g[0].decode() if isinstance(g[0], bytes) else g[0]
                      for g in gene[:]]
        print(f"[stage 0b] {len(gene_names)} genes")

        # Stream in chunks to avoid loading 50M rows at once
        all_rows = []
        chunk = 5_000_000
        kept = 0
        for start in range(0, total, chunk):
            end = min(start + chunk, total)
            ex = exp[start:end]
            x = ex["x"].astype(np.int32)
            y = ex["y"].astype(np.int32)
            count = ex["count"].astype(np.int32)
            g_idx = ex["geneID"].astype(np.int32)

            # Pull cell label from mask
            in_bounds = (x >= 0) & (x < W) & (y >= 0) & (y < H)
            cell_id = np.zeros(len(x), dtype=np.uint32)
            cell_id[in_bounds] = mask[y[in_bounds], x[in_bounds]]
            overlaps_nucleus = (cell_id > 0).astype(np.uint8)

            # Expand by count (each row in .gef is a (x,y,gene,count); proseg wants 1 row per molecule)
            # Use np.repeat for vectorised expansion
            n_mol = int(count.sum())
            x_e = np.repeat(x, count)
            y_e = np.repeat(y, count)
            g_e = np.repeat(g_idx, count)
            c_e = np.repeat(cell_id, count)
            o_e = np.repeat(overlaps_nucleus, count)

            df = pd.DataFrame({
                "transcript_id": np.arange(kept, kept + n_mol, dtype=np.uint64),
                "x_location": x_e.astype(np.float32),
                "y_location": y_e.astype(np.float32),
                "z_location": np.zeros(n_mol, dtype=np.float32),
                "gene": [gene_names[i] for i in g_e],
                "cell_id": c_e.astype(np.uint32),
                "overlaps_nucleus": o_e.astype(np.uint8),
            })
            all_rows.append(df)
            kept += n_mol
            print(f"[stage 0b]   chunk {start:,}-{end:,}  ({kept:,} molecules so far)")

    full = pd.concat(all_rows, ignore_index=True)
    print(f"[stage 0b] total molecules: {len(full):,}")
    print(f"[stage 0b]   assigned to cell : {(full['cell_id']>0).sum():,} ({100*(full['cell_id']>0).mean():.1f}%)")
    print(f"[stage 0b]   in nucleus       : {(full['overlaps_nucleus']==1).sum():,}")

    print(f"[stage 0b] writing {out_csv} ...")
    full.to_csv(out_csv, index=False, compression="gzip")
    print(f"[stage 0b] wrote {out_csv}  ({out_csv.stat().st_size/1e6:.1f} MB)")


def run_proseg(work_dir, proseg_binary, burnin=100, samples=100):
    """Run proseg with the validated settings. Output filenames are BARE (no path prefix)."""
    cmd = [
        proseg_binary,
        "transcripts.csv.gz",
        "--transcript-id-column", "transcript_id",
        "-x", "x_location", "-y", "y_location", "-z", "z_location",
        "--gene-column", "gene",
        "--cell-id-column", "cell_id",
        "--cell-id-unassigned", "0",
        "--compartment-column", "overlaps_nucleus",
        "--compartment-nuclear", "1",
        "--ignore-z-coord",
        "--burnin-samples", str(burnin),
        "--samples", str(samples),
        "--output-path", "output/",
        "--output-counts", "proseg-counts.mtx.gz",
        "--output-cell-metadata", "proseg-cell-metadata.csv.gz",
        "--output-gene-metadata", "proseg-gene-metadata.csv.gz",
        "--output-transcript-metadata", "proseg-transcript-metadata.csv.gz",
        "--output-cell-polygons", "proseg-cell-polygons.geojson.gz",
    ]
    print(f"[stage 0b] running proseg (working dir: {work_dir})")
    print(f"[stage 0b] cmd: {' '.join(cmd)}")
    rc = subprocess.call(cmd, cwd=work_dir)
    if rc != 0:
        sys.exit(f"proseg failed with exit code {rc}")


def load_proseg_to_h5ad(proseg_out_dir, out_h5ad):
    """Load proseg outputs (counts mtx + cell metadata) into an AnnData and save."""
    import anndata as ad
    from scipy.io import mmread

    counts = mmread(str(proseg_out_dir / "proseg-counts.mtx.gz")).tocsr()
    meta = pd.read_csv(proseg_out_dir / "proseg-cell-metadata.csv.gz")
    genes = pd.read_csv(proseg_out_dir / "proseg-gene-metadata.csv.gz")
    print(f"[stage 0b] loaded counts: {counts.shape}  cell meta: {meta.shape}  gene meta: {genes.shape}")

    a = ad.AnnData(X=counts, obs=meta, var=genes)
    a.obs_names = a.obs.get("cell", pd.Series(range(a.n_obs))).astype(str).values
    a.var_names = a.var.get("gene", a.var.iloc[:, 0]).astype(str).values

    # Build obsm['spatial'] from centroids
    if {"centroid_x", "centroid_y"}.issubset(a.obs.columns):
        a.obsm["spatial"] = a.obs[["centroid_x", "centroid_y"]].values.astype(np.float32)
    print(f"[stage 0b] AnnData: {a.shape}")
    a.write_h5ad(out_h5ad)
    print(f"[stage 0b] wrote {out_h5ad}  ({out_h5ad.stat().st_size/1e6:.1f} MB)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--skip-csv", action="store_true",
                        help="Skip rebuilding transcripts.csv.gz (if already exists).")
    parser.add_argument("--skip-proseg", action="store_true",
                        help="Skip the proseg call (if outputs already exist).")
    args = parser.parse_args()

    cfg = load_config(args.config)
    up = cfg.get("upstream", {})
    chip_id = up.get("chip_id", cfg["roi_id"])
    out_dir = Path(up.get("upstream_out_dir", f"outputs/{chip_id}_upstream"))
    out_dir.mkdir(parents=True, exist_ok=True)
    proseg_out = out_dir / "output"; proseg_out.mkdir(exist_ok=True)

    mask_path = out_dir / "nuclei_mask.npy.gz"
    csv_path  = out_dir / "transcripts.csv.gz"
    h5ad_path = out_dir / "proseg_full_raw.h5ad"

    proseg_bin = up.get("proseg_binary", "proseg")
    tissue_gef = Path(up.get("tissue_gef", ""))

    print(f"[stage 0b] chip_id         : {chip_id}")
    print(f"[stage 0b] mask            : {mask_path}")
    print(f"[stage 0b] tissue.gef      : {tissue_gef}")
    print(f"[stage 0b] proseg binary   : {proseg_bin}")
    print(f"[stage 0b] out h5ad        : {h5ad_path}")

    if not args.skip_csv:
        if not mask_path.exists():
            sys.exit(f"nuclei mask not found at {mask_path}. Run stage 0 (Cellpose) first.")
        if not tissue_gef.exists():
            sys.exit(f"tissue.gef not found at {tissue_gef}. Set upstream.tissue_gef in config.")
        print("[stage 0b] loading mask ...")
        with gzip.open(mask_path, "rb") as f:
            mask = np.load(f)
        build_transcripts_csv(tissue_gef, mask, csv_path)

    if not args.skip_proseg:
        if not csv_path.exists():
            sys.exit(f"transcripts.csv.gz not found at {csv_path}.")
        run_proseg(out_dir, proseg_bin,
                   burnin=up.get("proseg_burnin_samples", 100),
                   samples=up.get("proseg_samples", 100))

    load_proseg_to_h5ad(proseg_out, h5ad_path)

    print(f"\n[stage 0b] DONE. proseg_full_raw.h5ad written to {h5ad_path}")
    print(f"[stage 0b] Update your ROI YAML so input.proseg_h5ad points to this file,")
    print(f"[stage 0b] then run stages 1-9 normally:")
    print(f"[stage 0b]   python pipeline/run_all.py --config {args.config} --from-stage 1")


if __name__ == "__main__":
    main()
