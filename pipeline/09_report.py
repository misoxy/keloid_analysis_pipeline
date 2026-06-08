"""Stage 9: build a self-contained HTML report from all outputs.

Reads:  outputs/<roi_id>/figures/*.png and tables/*.csv
Writes: outputs/<roi_id>/report.html

Design:
  Images are embedded as base64 so the HTML is one portable file. The lab
  member can email or Dropbox it without losing the figures.

  Sections mirror the pipeline stages. Top-of-report summary lists key
  numbers: cell counts per type, IgG count, K=10 significant pair count.

Run:
    python pipeline/09_report.py --config configs/strip_01.yaml
"""
from pathlib import Path
import argparse
import base64
import warnings
warnings.filterwarnings("ignore")

import pandas as pd

import sys; sys.path.insert(0, str(Path(__file__).parent))
from io_utils import load_config, out_paths


def img_b64(path):
    if not path.exists():
        return f"<i>missing: {path.name}</i>"
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode("ascii")
    return f'<img src="data:image/png;base64,{data}" style="max-width:100%; border:1px solid #ddd; margin: 8px 0;">'


def csv_table(path, max_rows=20):
    if not path.exists():
        return f"<i>missing: {path.name}</i>"
    df = pd.read_csv(path)
    return df.head(max_rows).to_html(index=False, classes="t", border=0,
                                     float_format=lambda x: f"{x:.3f}")


CSS = """
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, sans-serif;
         max-width: 1100px; margin: 30px auto; color: #222; padding: 0 20px; }
  h1 { color: #c0392b; border-bottom: 2px solid #c0392b; padding-bottom: 8px; }
  h2 { color: #2c3e50; margin-top: 36px; border-bottom: 1px solid #ddd; padding-bottom: 4px; }
  h3 { color: #555; margin-top: 24px; }
  table.t { border-collapse: collapse; font-size: 13px; margin: 10px 0; }
  table.t th, table.t td { padding: 4px 10px; border-bottom: 1px solid #eee; text-align: left; }
  table.t th { background: #f5f5f5; }
  .summary { background: #fffbe5; padding: 14px; border-left: 4px solid #f39c12; }
  code { background: #f4f4f4; padding: 2px 4px; border-radius: 3px; }
</style>
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    paths = out_paths(cfg)
    title = cfg.get("report", {}).get("title", f"{cfg['roi_id']} report")

    # Pull a few headline numbers
    summary_html = ""
    sizes_csv = paths["tab"] / "05_igg_set_sizes.csv"
    if sizes_csv.exists():
        sizes = pd.read_csv(sizes_csv)
        summary_html += "<p><b>IgG candidate set sizes:</b></p>" + sizes.to_html(index=False, classes="t", border=0)
    sig_csv = paths["tab"] / "06_significant_pairs_K10.csv"
    if sig_csv.exists():
        n_sig = len(pd.read_csv(sig_csv))
        summary_html += f"<p><b>Neighbourhood significant pairs at K=10:</b> {n_sig}</p>"

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>{title}</title>{CSS}</head><body>

<h1>{title}</h1>
<p>ROI: <code>{cfg['roi_id']}</code> &nbsp; | &nbsp;
   coords: x[{cfg['strip_coords']['x_min']}, {cfg['strip_coords']['x_max']}],
           y[{cfg['strip_coords']['y_min']}, {cfg['strip_coords']['y_max']}]
</p>

<div class="summary"><b>Summary</b>{summary_html}</div>

<h2>Stage 3 -BANKSY lambda sweep</h2>
{img_b64(paths['fig'] / 'banksy_lambda_sweep.png')}

<h2>Stage 4 -Annotation</h2>
<h3>4a Broad cell types</h3>
{img_b64(paths['fig'] / '04a_broad_celltype_spatial.png')}
<h3>4b Fibroblast subtype proposal</h3>
{img_b64(paths['fig'] / '04b_fibroblast_subtype_spatial.png')}
<p>Marker panel coverage:</p>
{csv_table(paths['tab'] / '04b_fibroblast_marker_presence.csv', max_rows=10)}
<h3>4c Vessel marker presence</h3>
{csv_table(paths['tab'] / '04c_vessel_marker_presence.csv', max_rows=10)}
<h3>4d Detailed annotation</h3>
{img_b64(paths['fig'] / '04d_detailed_celltype_v1_spatial.png')}
<p>Per-label marker fractions (top 20 labels):</p>
{csv_table(paths['tab'] / '04d_detailed_label_marker_fraction.csv', max_rows=20)}

<h2>Stage 5 -IgG-rich candidate identification</h2>
{img_b64(paths['fig'] / '05_igg_spatial.png')}
<p>Sub-cluster IGHG1+ enrichment:</p>
{csv_table(paths['tab'] / '05_subcluster_ighg1_enrichment.csv', max_rows=15)}
<p>Marker fraction per candidate set:</p>
{csv_table(paths['tab'] / '05_igg_set_marker_fraction.csv', max_rows=20)}
<p>Sensitivity: alternative IgG definitions:</p>
{csv_table(paths['tab'] / '05_igg_sensitivity_definitions.csv', max_rows=10)}

<h2>Stage 6 -All-pair neighbourhood enrichment (K=10, spatial-bin null)</h2>
{img_b64(paths['fig'] / '06_neighbourhood_K10_heatmap.png')}
<p>IgG focal row (K=10):</p>
{csv_table(paths['tab'] / '06_IgG_focal_row_K10.csv', max_rows=15)}
<p>Significant pairs (top 25 by enrichment):</p>
{csv_table(paths['tab'] / '06_significant_pairs_K10.csv', max_rows=25)}

<h2>Stage 7 -Distance to vessel</h2>
{img_b64(paths['fig'] / '07_distance_boxplot.png')}
{csv_table(paths['tab'] / '07_distance_test.csv', max_rows=20)}

<h2>Stage 8 - Local stromal niche around IgG</h2>
{img_b64(paths['fig'] / '08_near_vs_far_spatial.png')}
<p>Niche program test (near vs far from IgG):</p>
{csv_table(paths['tab'] / '08_niche_program_test.csv', max_rows=15)}

<p style="color:#888; margin-top:40px; font-size:12px;">
This report contains only figures and statistical tables produced by the
pipeline. Interpretation is left to the reader. See
<code>provenance.json</code> in the output folder for the exact config and
git commit used to generate this report.
</p>
</body></html>"""

    paths["report"].write_text(html)
    print(f"[stage 9] wrote: {paths['report']}  ({paths['report'].stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
