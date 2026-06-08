"""Orchestrator: run all pipeline stages in order from one YAML config.

Supports re-running from any stage. Each stage is a separate python script
that reads/writes intermediate h5ads in outputs/<roi_id>/, so re-running
from stage N picks up where stage N-1 left off.

Usage:
    # Full pipeline
    python pipeline/run_all.py --config configs/strip_01.yaml

    # Re-run from stage 4 onward (e.g. you fixed annotation code)
    python pipeline/run_all.py --config configs/strip_01.yaml --from-stage 4

    # Run just one stage
    python pipeline/run_all.py --config configs/strip_01.yaml --from-stage 6 --to-stage 6

Stages:
    1 = subset strip
    2 = QC + normalise
    3 = BANKSY
    4 = annotate cells
    5 = IgG detection + sensitivity
    6 = all-pair neighbourhood
    7 = distance to vessel
    8 = local stromal niche
    9 = HTML report
"""
from pathlib import Path
import argparse
import subprocess
import sys
import time

THIS_DIR = Path(__file__).parent

STAGES = [
    # Upstream stages (slow, server-intended). Skipped by default in run_all
    # because most users will start from a pre-existing proseg h5ad. To run
    # them, use --from-stage 0.
    (0,  "00_cellpose_segment.py",   "Cellpose nuclei segmentation (server)"),
    (0.5, "00b_proseg_run.py",       "proseg cell segmentation (server)"),
    # Downstream stages (fast, per-ROI; what lab members typically run)
    (1, "01_subset_strip.py",       "subset strip"),
    (2, "02_qc_normalise.py",       "QC + normalise"),
    (3, "03_banksy_cluster.py",     "BANKSY lambda sweep"),
    (4, "04_annotate_cells.py",     "annotate cells"),
    (5, "05_igg_detection.py",      "IgG detection + sensitivity"),
    (6, "06_neighbourhood.py",      "all-pair neighbourhood"),
    (7, "07_distance_to_vessel.py", "distance to vessel"),
    (8, "08_local_niche.py",        "local stromal niche"),
    (9, "09_report.py",             "build HTML report"),
    (11, "11_tls_detect.py",        "tertiary lymphoid structure detection"),
]


def run_stage(n, script, label, config_path):
    print("\n" + "=" * 70)
    print(f"STAGE {n}: {label}  ({script})")
    print("=" * 70)
    t0 = time.time()
    rc = subprocess.call(
        [sys.executable, str(THIS_DIR / script), "--config", str(config_path)],
    )
    dt = time.time() - t0
    if rc != 0:
        print(f"\nSTAGE {n} FAILED with exit code {rc} after {dt:.1f}s.")
        print(f"Fix the issue and re-run:")
        print(f"  python pipeline/run_all.py --config {config_path} --from-stage {n}")
        sys.exit(rc)
    print(f"STAGE {n} done in {dt:.1f}s")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="path to ROI YAML config")
    parser.add_argument("--from-stage", type=float, default=1,
                        help="first stage to run. Default 1 (skip upstream Cellpose+proseg). "
                             "Use 0 to include Cellpose, 0.5 to include proseg.")
    parser.add_argument("--to-stage", type=float, default=11,
                        help="last stage to run (default 11; covers TLS detection after report)")
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    if not config_path.exists():
        sys.exit(f"config not found: {config_path}")

    print(f"config       : {config_path}")
    print(f"from stage   : {args.from_stage}")
    print(f"to stage     : {args.to_stage}")

    t_start = time.time()
    for n, script, label in STAGES:
        if n < args.from_stage or n > args.to_stage:
            continue
        run_stage(n, script, label, config_path)

    # Write provenance once at the end of a full run (or after the last stage)
    if args.to_stage == 9:
        # Lazy import to avoid hard dep when run_all is the only thing called
        sys.path.insert(0, str(THIS_DIR))
        from io_utils import load_config, out_paths, write_provenance
        cfg = load_config(config_path)
        paths = out_paths(cfg)
        write_provenance(cfg, config_path, paths,
                         extra={"total_runtime_seconds": round(time.time() - t_start, 1)})

    print(f"\n{'=' * 70}")
    print(f"PIPELINE COMPLETE in {time.time() - t_start:.1f}s")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
