# keloid_analysis_pipeline

Spatial transcriptomics analysis pipeline for keloid Stereo-seq data
(chip C04029G4, sample 5BTL_5).

The pipeline takes proseg-segmented cells inside a tissue strip, runs BANKSY
clustering at multiple lambda values, annotates clusters by marker genes,
and quantifies cell-cell spatial neighbourhoods. The output is an annotated
AnnData file plus an HTML report per strip.

The eventual biological question is which immune cell types sit next to
which fibroblast subtypes and which cells flank blood vessels, in keloid
tissue at single-cell spatial resolution.

## Repository layout

```
keloid_analysis_pipeline/
├── README.md              this file
├── environment.yml        conda environment spec
├── .gitignore             files git should never track
│
├── configs/
│   ├── strip_01.yaml      parameters for one strip (one YAML per ROI)
│   └── skin_markers.yaml  marker genes per skin cell type
│
├── pipeline/              the actual scripts
│   ├── 01_subset_strip.py     proseg AnnData -> one strip subset
│   ├── 02_qc_normalise.py     filter cells, HVG, normalise
│   ├── 03_banksy_cluster.py   run BANKSY at a sweep of lambda values
│   ├── 04_annotate.py         label clusters by marker gene expression
│   ├── 05_neighbourhood.py    compute cell-cell spatial interactions
│   ├── 06_report.py           build an HTML report
│   └── run_all.py             orchestrator that runs 01-06 from a config
│
├── notebooks/
│   ├── exploration/       dev notebooks (e.g. strip coordinate picker)
│   └── reports/           rendered HTML reports go here
│
└── outputs/               per-strip outputs (gitignored, never committed)
    └── <roi_id>/
```

## Setup

1. Clone the repo.
2. Build the conda env once:
   ```
   conda env create -f environment.yml
   conda activate keloid_pipeline
   ```
3. Make sure the raw data is on the machine (see paths in the configs).

## Running a strip

```
python pipeline/run_all.py --config configs/strip_01.yaml
```

To analyse a new strip, copy an existing config, change `roi_id` and
`strip_coords`, and rerun the same command.

## Status

Currently in Phase 0 (validation on a single strip using the existing
proseg cells). Cellpose integration is added later, once the BANKSY
analysis is locked.
