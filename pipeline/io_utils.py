"""Shared IO helpers: config loading, marker panel loading, provenance writing.

Used by every stage script and by run_all.py. Keeps the per-stage scripts
focused on biology, not on YAML parsing or path bookkeeping.
"""
from pathlib import Path
import subprocess
import json
import datetime
import yaml


def load_config(path):
    """Read a ROI YAML config and return it as a dict."""
    with open(path) as f:
        return yaml.safe_load(f)


def load_marker_panels(path):
    """Read the marker panel library YAML and return the full dict.

    Caller picks the panel family by key, e.g. cfg['fibroblast_subtype_markers'].
    """
    with open(path) as f:
        return yaml.safe_load(f)


def patch_log1p_base(adata):
    """Workaround for scanpy 1.9 + anndata h5ad round-trip that drops the
    'base' key from uns['log1p']. Several scanpy tools (rank_genes_groups,
    score_genes) raise KeyError without it. Call right after read_h5ad.
    """
    if "log1p" in adata.uns and "base" not in adata.uns["log1p"]:
        adata.uns["log1p"]["base"] = None
    return adata


def out_paths(cfg):
    """Standard set of output paths for one ROI. All stages use these names."""
    out = Path(cfg["output_dir"])
    fig = out / "figures"
    tab = out / "tables"
    out.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    tab.mkdir(parents=True, exist_ok=True)
    return {
        "out":   out,
        "fig":   fig,
        "tab":   tab,
        "strip":        out / "adata_strip.h5ad",
        "normalised":   out / "adata_normalised.h5ad",
        "banksy":       out / "adata_banksy.h5ad",
        "annotated":    out / "adata_annotated.h5ad",
        "with_igg":     out / "adata_with_igg.h5ad",
        "atlas":        out / "adata_atlas_full.h5ad",
        "provenance":   out / "provenance.json",
        "report":       out / "report.html",
    }


def git_commit_hash():
    """Best-effort git commit hash for the repo. Returns 'unknown' if not in a git repo."""
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
            cwd=Path(__file__).parent,
        )
        return out.decode().strip()
    except Exception:
        return "unknown"


def write_provenance(cfg, config_path, paths, extra=None):
    """Write provenance.json with config + git hash + timestamp + extras.

    Call this at the end of run_all.py so any figure or table in the output
    can be traced back to the exact code state and config that produced it.
    """
    rec = {
        "roi_id":        cfg.get("roi_id"),
        "config_path":   str(config_path),
        "config":        cfg,
        "git_commit":    git_commit_hash(),
        "timestamp_utc": datetime.datetime.utcnow().isoformat() + "Z",
    }
    if extra:
        rec.update(extra)
    with open(paths["provenance"], "w") as f:
        json.dump(rec, f, indent=2, default=str)
    print(f"  provenance -> {paths['provenance']}")
