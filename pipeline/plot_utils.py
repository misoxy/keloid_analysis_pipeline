"""Shared plotting helpers used by multiple stages.

Keeps figure styling consistent across the pipeline so the lab member sees
one visual language across all stage outputs.
"""
import numpy as np
import matplotlib.pyplot as plt


def spatial_scatter(adata, label_key, title=None, filename=None,
                    size=3, alpha=0.85, highlight=None, palette=None):
    """Spatial scatter coloured by a categorical obs column.

    highlight: list of label names to draw larger and on top (e.g. ['IgG_producing']).
    palette: optional dict label -> hex colour; otherwise tab20.
    """
    xy = np.asarray(adata.obsm["spatial"])
    labels = adata.obs[label_key].astype(str)
    cats = sorted(labels.unique())
    if palette is None:
        cmap = plt.cm.tab20(np.linspace(0, 1, max(len(cats), 1)))
        palette = {c: cmap[i] for i, c in enumerate(cats)}

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.scatter(xy[:, 0], xy[:, 1], s=1, c="lightgrey", alpha=0.2, edgecolors="none")
    for lab in cats:
        if highlight and lab in highlight:
            continue
        m = (labels == lab).values
        ax.scatter(xy[m, 0], xy[m, 1], s=size, color=palette[lab],
                   alpha=alpha, edgecolors="none",
                   label=f"{lab} (n={int(m.sum()):,})")
    if highlight:
        for lab in highlight:
            m = (labels == lab).values
            if m.sum() == 0:
                continue
            ax.scatter(xy[m, 0], xy[m, 1], s=size * 5, c="red",
                       edgecolors="black", linewidths=0.3, alpha=0.95,
                       label=f"{lab} (n={int(m.sum())})")
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title or label_key)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=8, markerscale=3, frameon=False)
    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=180, bbox_inches="tight")
    plt.close(fig)


def neighbour_heatmap(log2_df, qval_df, title="", filename=None, sig_q=0.05):
    """Heatmap of log2(observed/expected) with q<sig_q cells marked with *."""
    fig, ax = plt.subplots(
        figsize=(0.65 * log2_df.shape[1] + 4, 0.65 * log2_df.shape[0] + 3),
    )
    vmax = float(np.nanpercentile(np.abs(log2_df.values), 95)) or 1.0
    im = ax.imshow(log2_df.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(log2_df.shape[1]))
    ax.set_yticks(range(log2_df.shape[0]))
    ax.set_xticklabels(log2_df.columns, rotation=90)
    ax.set_yticklabels(log2_df.index)
    plt.colorbar(im, ax=ax, label="log2(observed / expected)")
    for i, focal in enumerate(log2_df.index):
        for j, neighbour in enumerate(log2_df.columns):
            if qval_df.loc[focal, neighbour] < sig_q:
                ax.text(j, i, "*", ha="center", va="center",
                        color="black", fontsize=9)
    ax.set_xlabel("neighbour cell type")
    ax.set_ylabel("focal cell type")
    ax.set_title(title)
    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=180, bbox_inches="tight")
    plt.close(fig)


def distance_boxplot(distances, labels, title="", filename=None, vline=None):
    """Per-label boxplot of distances, sorted by median (closest at top)."""
    import pandas as pd
    df = pd.DataFrame({"label": labels, "d": distances}).dropna()
    order = df.groupby("label")["d"].median().sort_values().index.tolist()
    fig, ax = plt.subplots(figsize=(11, max(4, 0.4 * len(order))))
    data = [df.loc[df["label"] == lab, "d"].values for lab in order]
    bp = ax.boxplot(data, vert=False, positions=np.arange(len(order)),
                    widths=0.6, patch_artist=True, showfliers=False)
    for p in bp["boxes"]:
        p.set_facecolor("#9ecae1"); p.set_alpha(0.7)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(order)
    ax.set_xlabel("distance (px)")
    ax.set_title(title)
    if vline is not None:
        ax.axvline(vline, color="red", ls="--", label="overall median")
        ax.legend()
    plt.tight_layout()
    if filename:
        plt.savefig(filename, dpi=180, bbox_inches="tight")
    plt.close(fig)
