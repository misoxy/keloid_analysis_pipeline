# keloid_analysis_pipeline

Spatial transcriptomics analysis pipeline for keloid Stereo-seq data, chip C04029G4. Takes proseg-segmented cells, subsets to a tissue strip, runs BANKSY clustering, annotates cells by marker panels, and quantifies which cell types sit next to which.

The biological goal is to map immune cells (especially IgG-rich antibody-secreting candidates) against fibroblast subtypes and vessels in keloid tissue at single-cell spatial resolution, something bin100 deconvolution cannot reliably do because rare immune signals are masked by abundant fibroblast collagen transcripts.

## What you get out

Per ROI, in `outputs/<roi_id>/`:

- `adata_atlas_full.h5ad` — the final annotated AnnData with all label columns.
- `figures/` — one PNG per analysis stage, numbered by stage.
- `tables/` — CSVs for every statistical test (neighbour enrichment, distance tests, niche analysis).
- `report.html` — single-file self-contained HTML report you can read on Dropbox without re-running anything.
- `provenance.json` — which config, which git commit, which input file, when it ran. So you can trace any figure back to the data.

## Quickstart for a lab member

You need three things: a clone of this repo, conda installed, and a path to a proseg AnnData (`.h5ad`) for the chip.

```bash
# 1. Clone the repo
git clone https://github.com/misoxy/keloid_analysis_pipeline.git
cd keloid_analysis_pipeline

# 2. Build the conda env (one-time, ~10 min)
conda env create -f environment.yml
conda activate keloid_pipeline

# 3. Edit configs/strip_01.yaml so the proseg_h5ad path points to your file,
#    then run the pipeline end-to-end
python pipeline/run_all.py --config configs/strip_01.yaml
```

That writes everything into `outputs/strip_01/`. Open `outputs/strip_01/report.html` in a browser.

## Running a new ROI

```bash
# 1. Copy the template
cp configs/template_roi.yaml configs/strip_02.yaml

# 2. Edit strip_02.yaml:
#    - change roi_id to "strip_02"
#    - change strip_coords (x_min, x_max, y_min, y_max in chip pixel space)
#    - change output_dir to outputs/strip_02
#    - leave QC / BANKSY / panel paths alone unless you have a reason

# 3. Run
python pipeline/run_all.py --config configs/strip_02.yaml
```

## Re-running a single stage

If BANKSY succeeded but annotation failed, fix the annotation code, then:

```bash
# Re-run only from stage 4 onward
python pipeline/run_all.py --config configs/strip_01.yaml --from-stage 4

# Or just one stage
python pipeline/run_all.py --config configs/strip_01.yaml --from-stage 6 --to-stage 6

# Or directly
python pipeline/06_neighbourhood.py --config configs/strip_01.yaml
```

All stages read intermediate `.h5ad` files from `outputs/<roi_id>/`, so you can resume anywhere.

## Pipeline stages

| # | Script | Reads | Writes |
|---|---|---|---|
| 1 | `01_subset_strip.py` | proseg full chip h5ad | `adata_strip.h5ad` |
| 2 | `02_qc_normalise.py` | strip h5ad | `adata_normalised.h5ad` |
| 3 | `03_banksy_cluster.py` | normalised h5ad | `adata_banksy.h5ad`, `figures/03_banksy_lambda_sweep.png` |
| 4 | `04_annotate_cells.py` | banksy h5ad + marker panels | `adata_annotated.h5ad`, marker tables, spatial PNGs |
| 5 | `05_igg_detection.py` | annotated h5ad | `adata_with_igg.h5ad`, IgG sensitivity table |
| 6 | `06_neighbourhood.py` | igg h5ad | neighbourhood CSVs, K=10 heatmap PNG |
| 7 | `07_distance_to_vessel.py` | igg h5ad | distance CSVs, distance boxplot PNG |
| 8 | `08_local_niche.py` | igg h5ad | niche test CSVs, near-vs-far PNGs |
| 9 | `09_report.py` | all of the above | `report.html` (self-contained, base64 images) |

## Configuration

Each ROI is described by one YAML file in `configs/`. The template
`configs/template_roi.yaml` documents every field with inline comments.

Marker panels are in `configs/marker_panels.yaml`, organised into four panel families:

- `broad_markers` — basal keratinocyte, suprabasal keratinocyte, pan-fibroblast, myofibroblast, endothelial, lymphatic endothelial, pericyte/smooth muscle, macrophage, mast cell, T cell, B cell.
- `fibroblast_subtype_markers` — mesenchymal / secretory papillary / secretory reticular / pro-inflammatory / myofibroblast.
- `vessel_markers` — endothelial / lymphatic endothelial / pericyte-smooth muscle.
- `immune_markers` — macrophage / mast cell / T cell / B cell / IgG_producing.

To add a new cell type, add a key under the relevant panel family. Genes are HGNC symbols matching the AnnData `var_names`.

## Repository layout

```
keloid_analysis_pipeline/
├── README.md                 this file
├── environment.yml           conda env spec (named keloid_pipeline)
├── .gitignore                outputs/, *.h5ad, etc.
│
├── configs/
│   ├── template_roi.yaml     copy + edit for a new ROI
│   ├── strip_01.yaml         working example
│   └── marker_panels.yaml    all marker panels in one file
│
├── pipeline/
│   ├── __init__.py
│   ├── run_all.py            orchestrator (--from-stage / --to-stage)
│   ├── io_utils.py           config loading, h5ad IO, provenance
│   ├── plot_utils.py         shared plotting helpers
│   ├── stats_utils.py        BH-FDR, neighbour enrichment, spatial-bin perm
│   ├── 01_subset_strip.py
│   ├── 02_qc_normalise.py
│   ├── 03_banksy_cluster.py
│   ├── 04_annotate_cells.py
│   ├── 05_igg_detection.py
│   ├── 06_neighbourhood.py
│   ├── 07_distance_to_vessel.py
│   ├── 08_local_niche.py
│   └── 09_report.py
│
├── notebooks/
│   └── exploration/          dev notebooks
│
└── outputs/                  per-ROI outputs (gitignored)
    └── <roi_id>/
        ├── adata_strip.h5ad
        ├── adata_normalised.h5ad
        ├── adata_banksy.h5ad
        ├── adata_annotated.h5ad
        ├── adata_with_igg.h5ad
        ├── adata_atlas_full.h5ad
        ├── figures/
        ├── tables/
        ├── provenance.json
        └── report.html
```

## Expected runtimes (Apple Silicon laptop, ~8K cells per strip)

| Stage | Time |
|---|---|
| 1 subset | seconds |
| 2 QC/normalise | ~30 sec |
| 3 BANKSY | ~5-10 min |
| 4 annotate | ~30 sec |
| 5 IgG detection + sensitivity | ~30 sec |
| 6 neighbourhood (K=5,10,20, 500 perms each, spatial-bin null) | ~3-5 min |
| 7 distance-to-vessel | ~30 sec |
| 8 local niche | ~30 sec |
| 9 report | seconds |
| **Total per ROI** | **~10-20 min** |

Cellpose + proseg on the full chip (the upstream segmentation) is hours on a laptop, minutes on the server with GPU. They live outside this pipeline for now and produce the input `proseg_h5ad`.

## Caveats

- Marker-score-based subtyping is a *proposal*, not a classification. Cells whose top panel score is below the z threshold are flagged `*_uncertain` instead of force-classified.
- `IgG_rich_candidate` is the honest label for the focal IgG aggregate. Stereo-seq mRNA cannot resolve mature plasma cells from antibody-internalising macrophages without protein validation (CD138, MZB1, IgG, CD45/CD79A on IF/IHC).
- A single strip is one observation. Reproducibility requires running this pipeline on multiple strips and comparing.
- Permutation tests use spatial-bin null by default (preserves coarse tissue geography), not global label shuffle (too liberal in spatially structured tissue).

## Troubleshooting

- **`No such file: proseg_h5ad`** — edit the `input.proseg_h5ad` line in your YAML config.
- **`AssertionError: cluster N not in LABEL_MAP`** — your BANKSY clusters were numbered differently than the locked map in `04_annotate_cells.py`. Look at the marker table printed by stage 4 and update the LABEL_MAP at the top of that script.
- **`KeyError: 'banksy_l0.2'`** — re-run stage 3 with the correct lambda values.
- **`ModuleNotFoundError: banksy`** — env not active. Run `conda activate keloid_pipeline`.
- **Stage 6 takes forever** — reduce `n_perm` in the config from 500 to 200 for a quick check.
