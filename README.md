# keloid_analysis_pipeline

This repository contains a reproducible pipeline for analysing keloid Stereo-seq data at single-cell spatial resolution. It starts from proseg-segmented cells, crops user-defined tissue strips, runs QC and BANKSY clustering, annotates major cell types using marker panels, and tests which cell populations are physically close to each other.

## Starting guide

### What this is for

Keloid scars contain abnormal interactions between fibroblasts and immune cells that drive scar overgrowth. Bulk RNA and coarse bin-level data cannot resolve these interactions because the abundant fibroblast collagen signal drowns out the rarer immune signal. This pipeline operates on proseg-segmented single cells with BANKSY clustering and spatial neighbourhood statistics, so you can ask quantitative questions such as "do IgG-producing cells sit next to a specific fibroblast subtype?" with permutation-based p-values and FDR control.

### Guide

One-time setup on your machine:

```bash
git clone https://github.com/misoxy/keloid_analysis_pipeline.git
cd keloid_analysis_pipeline
conda env create -f environment.yml
conda activate keloid_pipeline
```

To analyse a tissue strip:

1. Open `configs/strip_01.yaml` and change the `input.proseg_h5ad` line so it points at your proseg AnnData file. If you are picking a new strip on the chip, also change the `strip_coords` block.
2. Run:
   ```bash
   python pipeline/run_all.py --config configs/strip_01.yaml
   ```
3. When it finishes, open `outputs/strip_01/report.html` in a browser. That is your full result.

To analyse a new strip without overwriting strip_01, copy the config first:

```bash
cp configs/strip_01.yaml configs/strip_02.yaml
```

Then edit `strip_02.yaml` (change `roi_id`, `strip_coords`, `output_dir`), and run with `--config configs/strip_02.yaml`.

## Required input format

The pipeline reads an AnnData `.h5ad` file at the path you put in `input.proseg_h5ad`. The file must satisfy:

- `adata.X`: a (cells × genes) expression matrix, sparse or dense, raw integer counts (not pre-normalised).
- `adata.obs`: one row per cell. `adata.obs_names` should be unique cell IDs.
- `adata.var`: one row per gene. `adata.var_names` MUST be HGNC gene symbols matching the marker panels in `configs/marker_panels.yaml`. If your genes are Ensembl IDs, convert to symbols first.
- `adata.obsm["spatial"]`: an (n_cells × 2) array of (x, y) centroid coordinates in chip pixel space.
- Source: typically produced by `proseg` on proseg-segmented cells. Other segmenters work as long as the four requirements above are met.

If your file is missing any of these, an early pipeline stage will fail with a clear error message naming the missing field.

## What you get out (output)

Per ROI, in `outputs/<roi_id>/`:

- `adata_atlas_full.h5ad`: the final annotated AnnData with every label column from stages 1-8. Load this in Python (`anndata.read_h5ad`) for any ad-hoc analysis beyond the report.
- `figures/`: one PNG per analysis stage, numbered by stage (e.g. `04a_broad_celltype_spatial.png`, `06_neighbourhood_K10_heatmap.png`).
- `tables/`: CSVs for every statistical test (neighbour enrichment, distance tests, niche analysis, IgG sensitivity).
- `report.html`: single-file self-contained HTML report you can open in a browser without re-running anything. All figures embedded as base64, all tables rendered. Portable: send it via Dropbox / email and the recipient sees everything.
- `provenance.json`: which config, which git commit, which input file, when it ran, total runtime. Lets you trace any figure back to the data.

Intermediate AnnData files (`adata_strip.h5ad`, `adata_normalised.h5ad`, `adata_banksy.h5ad`, `adata_annotated.h5ad`, `adata_with_igg.h5ad`) are written between stages so you can resume with `--from-stage N` if a later stage fails.

## Reading the report

The HTML report sections correspond directly to pipeline stages. The numbers and figures are produced by the code; interpretation is left to the reader.

| Report section | What to look at |
|---|---|
| Summary | IgG candidate set sizes (Set D strict, Set A loose, etc.) and total significant neighbour pairs at K=10. Quick health-check. |
| Stage 3 BANKSY lambda sweep | Did the lambda sweep produce sensible cluster counts? Each lambda should give 3-8 clusters. |
| Stage 4 Annotation | Broad cell counts. The marker presence tables tell you which marker panels had missing genes in your dataset (low coverage = weak calls). The detailed marker fraction table is the per-label sanity check. |
| Stage 5 IgG identification | The selected sub-cluster's IGHG1+ fraction should be ≥ 50%. The set sizes show how many cells each definition gives. The marker fraction per set shows whether heavy + light chains are co-enriched. |
| Stage 6 Neighbourhood enrichment | Heatmap with `*` marking q < 0.05 pairs. The IgG focal row table tells you what sits next to IgG cells. |
| Stage 7 Distance to vessel | Boxplot of per-label distances. The test table tells you which labels live closer to vessels than chance. |
| Stage 8 Local stromal niche | Program-level table comparing near-IgG vs far-from-IgG stromal cells. Significant programs (q < 0.05) name candidate niche signatures. |

## Pipeline stages

| # | Script | Reads | Writes |
|---|---|---|---|
| 0 | `00_cellpose_segment.py` | H&E TIFF | `nuclei_mask.npy.gz` (slow, server-intended) |
| 0b | `00b_proseg_run.py` | mask + tissue.gef | `proseg_full_raw.h5ad` (slow, server-intended) |
| 1 | `01_subset_strip.py` | proseg full chip h5ad | `adata_strip.h5ad` |
| 2 | `02_qc_normalise.py` | strip h5ad | `adata_normalised.h5ad` |
| 3 | `03_banksy_cluster.py` | normalised h5ad | `adata_banksy.h5ad`, `figures/03_banksy_lambda_sweep.png` |
| 4 | `04_annotate_cells.py` | banksy h5ad + marker panels | `adata_annotated.h5ad`, marker tables, spatial PNGs |
| 5 | `05_igg_detection.py` | annotated h5ad | `adata_with_igg.h5ad`, IgG sensitivity table |
| 6 | `06_neighbourhood.py` | igg h5ad | neighbourhood CSVs, K=10 heatmap PNG |
| 7 | `07_distance_to_vessel.py` | igg h5ad | distance CSVs, distance boxplot PNG |
| 8 | `08_local_niche.py` | igg h5ad | niche test CSVs, near-vs-far PNGs |
| 9 | `09_report.py` | all of the above | `report.html` (self-contained, base64 images) |

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

## Configuration

Each ROI is described by one YAML file in `configs/`. The template `configs/template_roi.yaml` documents every field with inline comments.

Marker panels are in `configs/marker_panels.yaml`, organised into four panel families:

- `broad_markers`: basal keratinocyte, suprabasal keratinocyte, pan-fibroblast, myofibroblast, endothelial, lymphatic endothelial, pericyte/smooth muscle, macrophage, mast cell, T cell, B cell.
- `fibroblast_subtype_markers`: mesenchymal / secretory papillary / secretory reticular / pro-inflammatory / myofibroblast.
- `vessel_markers`: endothelial / lymphatic endothelial / pericyte-smooth muscle.
- `immune_markers`: macrophage / mast cell / T cell / B cell / IgG_producing.

To add a new cell type, add a key under the relevant panel family. Genes are HGNC symbols matching the AnnData `var_names`.

## Caveats

- Marker-score-based subtyping is a proposal, not a classification. Cells whose top panel score is below the z threshold are flagged `*_uncertain` rather than force-classified.
- `IgG_rich_candidate` is the honest label for the focal IgG aggregate. Stereo-seq mRNA cannot resolve mature plasma cells from antibody-internalising macrophages without protein validation (CD138, MZB1, IgG, CD45/CD79A on IF/IHC).
- A single strip is one observation. Reproducibility requires running this pipeline on multiple strips and comparing.
- Permutation tests use a spatial-bin null by default, which preserves coarse tissue geography, rather than a global label shuffle which would be too liberal in spatially structured tissue.

## Troubleshooting

- `No such file: proseg_h5ad`: edit the `input.proseg_h5ad` line in your YAML config.
- `AssertionError: cluster N not in LABEL_MAP`: your BANKSY clusters were numbered differently than the locked map. Look at the marker table printed by stage 4 and update the `broad_label_map` in your YAML.
- `KeyError: 'banksy_l0.2'`: re-run stage 3 with the correct lambda values.
- `ModuleNotFoundError: banksy`: the conda env is not active. Run `conda activate keloid_pipeline`.
- Stage 6 takes too long: reduce `n_perm` in the config from 500 to 200 for a quick check.
