"""
03_cluster_analysis.py
Hierarchical clustering, PCA, GABAergic co-expression, and cluster validation.

Reads batch-corrected, z-scored expression data from 02_preprocess.py and
produces five figures:

  fig1_donor_dendrogram.png     Ward's linkage dendrogram of donors, colored
                                by diagnosis. CPCC, silhouette, and ARI vs
                                diagnosis reported in title.

  fig2_pca.png                  PCA of donor expression profiles. Tests whether
                                the 48-gene panel separates SCZ from Control
                                in PC space.

  fig3_gabaergic_coexpression.png
                                Side-by-side Pearson correlation heatmaps of
                                GABAergic interneuron markers, SCZ vs Control.
                                Primary figure for Research Question 1: do
                                PV+, SST+, and VIP+ subtypes co-vary together
                                or independently?

  fig4_gene_heatmap.png         Expression heatmap with donors sorted by
                                diagnosis and genes clustered by Ward's linkage.

  fig5_cluster_validation.png   Silhouette score and ARI vs diagnosis across
                                k = 2 to 5 clusters.

Cluster validation metrics:
  CPCC (cophenetic correlation coefficient): measures how well the dendrogram
    preserves pairwise distances. >0.75 indicates good hierarchical structure.

  Silhouette score: measures within-cluster cohesion vs between-cluster
    separation. Range [-1, 1]; higher is better.

  ARI (Adjusted Rand Index): measures agreement between cluster assignments
    and true diagnosis labels. 1.0 = perfect, 0 = random. ARI near 0 in this
    analysis is consistent with published postmortem brain transcriptomics
    literature, where technical and biological covariates typically dominate
    expression variance.

References:
Guillozet-Bongaarts AL, et al. (2014). Altered gene expression in the
  dorsolateral prefrontal cortex of individuals with schizophrenia.
  Molecular Psychiatry, 19(4):478-485. doi:10.1038/mp.2013.30. PMID: 23528911.

Lanz TA, et al. (2019). Postmortem transcriptional profiling reveals widespread
  increase in inflammation in schizophrenia. Translational Psychiatry, 9(1):151.
  doi:10.1038/s41398-019-0494-6. PMID: 31123247. GEO: GSE53987.

Tan PN, Steinbach M, Kumar V (2006). Introduction to Data Mining.
  Chapter 8: Cluster Analysis. Pearson.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, cophenet, leaves_list
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, silhouette_score
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

DATA_DIR = Path(__file__).parent.parent / "data"
FIG_DIR = Path(__file__).parent.parent / "figures"
FIG_DIR.mkdir(exist_ok=True)

SCZ_COLOR = "#C1392B"
CTRL_COLOR = "#2980B9"

SUBTYPE_COLORS = {
    "PV+": "#7B2D8B",
    "SST+": "#2471A3",
    "VIP+": "#1E8449",
    "CB+": "#D4AC0D",
    "Pan-GABA": "#BA4A00",
    "Excitatory": "#7D3C98",
    "Oligodendrocyte": "#566573",
    "Risk/Signaling": "#909497",
    "Synaptic": "#B2BABB",
    "SCZ candidate": "#D5D8DC",
}

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    0.8,
    "font.family":       "DejaVu Sans",
    "savefig.dpi":       180,
    "savefig.bbox":      "tight",
    "savefig.facecolor": "white",
})


def load():
    expr = pd.read_csv(DATA_DIR / "expression_matrix.csv", index_col=0)
    meta = pd.read_csv(DATA_DIR / "donor_metadata.csv", index_col=0)
    genes = pd.read_csv(DATA_DIR / "gene_annotations.csv")
    meta = meta.reindex(expr.index)
    print(f"Loaded: {expr.shape}  (donors x genes)")
    return expr, meta, genes


# Fig 1: Dendrogram
def plot_dendrogram(expr, meta):
    """
    Ward's linkage hierarchical clustering of donors.
    Color bar below dendrogram shows diagnosis for each leaf.
    """
    print("  Fig 1: Dendrogram...")
    dist_vec = pdist(expr.values, metric="euclidean")
    Z = linkage(dist_vec, method="ward")
    cpcc, _ = cophenet(Z, dist_vec)
    labels_2 = fcluster(Z, t=2, criterion="maxclust")
    sil = silhouette_score(expr.values, labels_2)

    ari = None
    if "diagnosis" in meta.columns:
        diag = meta["diagnosis"].reindex(expr.index)
        true = (diag == "Schizophrenia").astype(int)
        ari = adjusted_rand_score(true.values, labels_2)

    # Two-panel: dendrogram on top, thin diagnosis color bar below
    fig = plt.figure(figsize=(14, 5.5))
    gs = gridspec.GridSpec(2, 1, height_ratios=[11, 1], hspace=0.04)
    ax_tree = fig.add_subplot(gs[0])
    ax_bar = fig.add_subplot(gs[1])

    dn = dendrogram(
        Z, ax=ax_tree,
        color_threshold=0,
        above_threshold_color="#AAAAAA",
        no_labels=True,
        link_color_func=lambda k: "#AAAAAA",
    )

    if "diagnosis" in meta.columns:
        leaf_order = dn["leaves"]
        colors = [
            SCZ_COLOR if diag.iloc[i] == "Schizophrenia" else CTRL_COLOR
            for i in leaf_order
        ]
        for j, col in enumerate(colors):
            ax_bar.add_patch(plt.Rectangle((j, 0), 1, 1, color=col, linewidth=0))
        ax_bar.set_xlim(0, len(colors))
        ax_bar.set_ylim(0, 1)
        ax_bar.axis("off")

        ax_bar.legend(
            handles=[
                mpatches.Patch(color=SCZ_COLOR, label="Schizophrenia"),
                mpatches.Patch(color=CTRL_COLOR, label="Control"),
            ],
            loc="lower center", ncol=2, bbox_to_anchor=(0.5, -1.8),
            frameon=True, edgecolor="#CCCCCC", fontsize=9
        )

    ari_str = f"  |  ARI vs diagnosis = {ari:.3f}" if ari is not None else ""
    ax_tree.set_title(
        f"Hierarchical Clustering of Donors - Ward's Linkage, Euclidean Distance\n"
        f"CPCC = {cpcc:.3f}  |  Silhouette (k=2) = {sil:.3f}{ari_str}",
        fontsize=11, pad=8
    )
    ax_tree.set_ylabel("Distance (Ward linkage)", fontsize=11)
    ax_tree.spines["top"].set_visible(False)
    ax_tree.spines["right"].set_visible(False)
    ax_tree.spines["bottom"].set_visible(False)
    ax_tree.tick_params(bottom=False)

    fig.savefig(FIG_DIR / "fig1_donor_dendrogram.png")
    plt.close()


# Fig 2: PCA
def plot_pca(expr, meta):
    """
    PCA of donor expression profiles.
    ARI from median split of PC1 tests whether PC1 separates diagnosis.
    """
    print("  Fig 2: PCA...")
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(expr.values)
    v1, v2 = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(7, 6))

    if "diagnosis" in meta.columns:
        diag = meta["diagnosis"].reindex(expr.index)
        for label, color in [("Schizophrenia", SCZ_COLOR), ("Control", CTRL_COLOR)]:
            mask = (diag == label).values
            ax.scatter(
                coords[mask, 0], coords[mask, 1],
                c=color, s=60, alpha=0.75,
                edgecolors="white", linewidths=0.5,
                label=f"{label} (n={mask.sum()})", zorder=3
            )
        # ARI from median split of PC1
        true = (diag == "Schizophrenia").astype(int).values
        pc1_labels = (coords[:, 0] > np.median(coords[:, 0])).astype(int)
        ari = adjusted_rand_score(true, pc1_labels)
        ax.text(0.03, 0.97, f"PC1 split ARI = {ari:.3f}",
                transform=ax.transAxes, fontsize=9, va="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor="#CCCCCC", alpha=0.9))

    ax.axhline(0, color="#DDDDDD", lw=0.8)
    ax.axvline(0, color="#DDDDDD", lw=0.8)
    ax.set_xlabel(f"PC1 - {v1:.1f}% of variance", fontsize=11)
    ax.set_ylabel(f"PC2 - {v2:.1f}% of variance", fontsize=11)
    ax.set_title(
        "Principal Component Analysis of Donor Expression Profiles\n"
        "Does the 48-gene expression profile separate diagnosis?",
        fontsize=11
    )
    ax.legend(frameon=True, edgecolor="#CCCCCC", fontsize=9, loc="upper right")

    fig.savefig(FIG_DIR / "fig2_pca.png")
    plt.close()


# Fig 3: GABAergic co-expression
def plot_gabaergic_coexpression(expr, genes, meta):
    """
    Side-by-side Pearson correlation heatmaps for GABAergic markers,
    SCZ vs Control. Primary figure addressing Research Question 1.

    Key finding: in controls, PVALB is relatively independent of SST+
    markers; in schizophrenia this independence is lost and all markers
    co-vary as a single module.
    """
    print("  Fig 3: GABAergic co-expression...")
    gaba = genes[genes["is_gabaergic"]]["gene"].tolist()
    gaba = [g for g in gaba if g in expr.columns]
    if len(gaba) < 3:
        print("  Not enough GABAergic markers - skipping")
        return

    gene_to_sub = dict(zip(genes["gene"], genes["subtype"]))

    fig, axes = plt.subplots(1, 2, figsize=(16, 6.5))
    fig.subplots_adjust(top=0.88, wspace=0.45)

    for ax, (diag, label, color) in zip(axes, [
        ("Schizophrenia", "Schizophrenia (n=48)", SCZ_COLOR),
        ("Control", "Control (n=55)", CTRL_COLOR),
    ]):
        if "diagnosis" in meta.columns:
            sub = expr[gaba].loc[meta["diagnosis"] == diag]
        else:
            sub = expr[gaba]

        if len(sub) < 3:
            ax.set_title(f"Insufficient data ({diag})")
            continue

        corr = sub.corr(method="pearson")
        order = leaves_list(linkage(pdist(corr.values), method="ward"))
        corr = corr.iloc[order, order]

        sns.heatmap(
            corr, ax=ax,
            cmap="RdBu_r", vmin=-1, vmax=1, center=0,
            annot=True, fmt=".2f",
            annot_kws={"size": 9},
            linewidths=0.4, linecolor="#EEEEEE",
            square=True,
            cbar_kws={"label": "Pearson r", "shrink": 0.75, "pad": 0.02},
            xticklabels=True, yticklabels=True,
        )

        for tick in ax.get_xticklabels():
            g = tick.get_text()
            tick.set_color(SUBTYPE_COLORS.get(gene_to_sub.get(g, ""), "#333"))
            tick.set_fontweight("bold")
            tick.set_fontsize(9)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

        for tick in ax.get_yticklabels():
            g = tick.get_text()
            tick.set_color(SUBTYPE_COLORS.get(gene_to_sub.get(g, ""), "#333"))
            tick.set_fontweight("bold")
            tick.set_fontsize(9)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

        ax.set_title(label, fontsize=12, color=color, fontweight="bold", pad=10)
        ax.set_xlabel("")
        ax.set_ylabel("")

    fig.suptitle(
        "GABAergic Interneuron Co-expression: Schizophrenia vs Control\n"
        "Do PV+, SST+, and VIP+ subtypes co-vary together or independently?",
        fontsize=12, y=0.97
    )

    patches = [
        mpatches.Patch(color=SUBTYPE_COLORS[k], label=k)
        for k in ["PV+", "SST+", "VIP+", "CB+", "Pan-GABA"]
    ]
    fig.legend(
        handles=patches,
        loc="lower center", ncol=5,
        title="Interneuron subtype", title_fontsize=9,
        fontsize=9, frameon=True, edgecolor="#CCCCCC",
        bbox_to_anchor=(0.5, -0.10)
    )

    fig.savefig(FIG_DIR / "fig3_gabaergic_coexpression.png", bbox_inches="tight")
    plt.close()


# Fig 4: Gene heatmap
def plot_gene_heatmap(expr, meta, genes):
    """
    Expression heatmap with donors sorted by diagnosis (Control then SCZ)
    and genes clustered by Ward's linkage. Gene labels colored by subtype.
    """
    print("  Fig 4: Gene heatmap...")
    gene_to_sub = dict(zip(genes["gene"], genes["subtype"]))

    if "diagnosis" in meta.columns:
        diag = meta["diagnosis"].reindex(expr.index)
        order = diag.sort_values().index
        expr_s = expr.loc[order]
        diag_s = diag.loc[order]
    else:
        expr_s = expr
        diag_s = None

    Z_genes = linkage(pdist(expr_s.T.values, metric="euclidean"), method="ward")
    gene_order = leaves_list(Z_genes)
    expr_plot = expr_s.T.iloc[gene_order]   # genes x donors

    n_genes = len(expr_plot)
    n_donors = len(expr_s)

    fig = plt.figure(figsize=(16, max(10, n_genes * 0.28)))
    gs = gridspec.GridSpec(
        2, 2,
        height_ratios=[0.5, n_genes],
        width_ratios=[n_donors, 1],
        hspace=0.02, wspace=0.02
    )
    ax_diag = fig.add_subplot(gs[0, 0])   # diagnosis bar
    ax_heat = fig.add_subplot(gs[1, 0])   # main heatmap
    ax_sub = fig.add_subplot(gs[1, 1])    # subtype color strip

    # Diagnosis bar
    if diag_s is not None:
        for j, d in enumerate(diag_s.values):
            col = SCZ_COLOR if d == "Schizophrenia" else CTRL_COLOR
            ax_diag.add_patch(plt.Rectangle((j - 0.5, 0), 1, 1, color=col, linewidth=0))
        n_ctrl = (diag_s == "Control").sum()
        n_scz = (diag_s == "Schizophrenia").sum()
        ax_diag.text(n_ctrl / 2, 0.5, "Control",
                     ha="center", va="center", fontsize=10,
                     color="white", fontweight="bold")
        ax_diag.text(n_ctrl + n_scz / 2, 0.5, "Schizophrenia",
                     ha="center", va="center", fontsize=10,
                     color="white", fontweight="bold")
    ax_diag.set_xlim(-0.5, n_donors - 0.5)
    ax_diag.set_ylim(0, 1)
    ax_diag.axis("off")
    ax_diag.set_title(
        "Gene Expression Heatmap - Batch-Corrected, Z-Scored\n"
        "Donors sorted by diagnosis  |  Genes clustered by Ward's linkage",
        fontsize=11, pad=6
    )

    # Main heatmap
    im = ax_heat.imshow(
        expr_plot.values,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-2.5, vmax=2.5,
        interpolation="nearest",
    )

    ax_heat.set_yticks(range(n_genes))
    ax_heat.set_yticklabels(expr_plot.index, fontsize=8)
    for tick, gene in zip(ax_heat.get_yticklabels(), expr_plot.index):
        tick.set_color(SUBTYPE_COLORS.get(gene_to_sub.get(gene, ""), "#333333"))

    ax_heat.set_xticks([])
    ax_heat.set_xlabel("Donors  (<- Control  |  Schizophrenia ->)", fontsize=10)

    if diag_s is not None:
        n_ctrl = (diag_s == "Control").sum()
        ax_heat.axvline(n_ctrl - 0.5, color="black", lw=1.5, ls="--", alpha=0.7)

    fig.canvas.draw()
    pos = ax_heat.get_position()
    cbar_ax = fig.add_axes([pos.x1 + 0.072, pos.y0, 0.015, pos.height])
    cb = fig.colorbar(im, cax=cbar_ax)
    cb.set_ticks([-2, -1, 0, 1, 2])
    cb.set_label("")
    cb.ax.set_xlabel("z-score", fontsize=9, labelpad=6)
    cb.ax.xaxis.set_label_position("bottom")

    # Subtype color strip (right side)
    for j, gene in enumerate(expr_plot.index):
        col = SUBTYPE_COLORS.get(gene_to_sub.get(gene, ""), "#D5D8DC")
        ax_sub.add_patch(plt.Rectangle((0, j), 1, 1, color=col, linewidth=0))
    ax_sub.set_xlim(0, 1)
    ax_sub.set_ylim(0, n_genes)
    ax_sub.axis("off")

    present = set(gene_to_sub.get(g, "SCZ candidate") for g in expr_plot.index)
    patches = [mpatches.Patch(color=v, label=k)
               for k, v in SUBTYPE_COLORS.items() if k in present]
    fig.legend(
        handles=patches,
        loc="upper right", bbox_to_anchor=(1.12, 0.92),
        title="Gene subtype", title_fontsize=8,
        fontsize=8, frameon=True, edgecolor="#CCCCCC"
    )

    fig.savefig(FIG_DIR / "fig4_gene_heatmap.png", bbox_inches="tight")
    plt.close()


# Fig 5: Cluster validation
def plot_cluster_validation(expr, meta):
    """
    Silhouette score (internal validity) and ARI vs diagnosis (external validity)
    for k = 2 to 5 clusters. Both metrics annotated on the plot.
    """
    print("  Fig 5: Cluster validation...")
    Z = linkage(pdist(expr.values, metric="euclidean"), method="ward")
    rows = []
    for k in [2, 3, 4, 5]:
        labels = fcluster(Z, t=k, criterion="maxclust")
        sil = silhouette_score(expr.values, labels)
        ari = None
        if "diagnosis" in meta.columns:
            true = (meta["diagnosis"].reindex(expr.index) == "Schizophrenia").astype(int)
            ari = adjusted_rand_score(true.values, labels)
        rows.append({"k": k, "Silhouette score": sil, "ARI vs diagnosis": ari})

    df = pd.DataFrame(rows)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle("Cluster Validation", fontsize=12, y=1.02)

    axes[0].plot(df["k"], df["Silhouette score"], "o-",
                 color="#2980B9", lw=2, ms=8, zorder=3)
    axes[0].fill_between(df["k"], df["Silhouette score"],
                         alpha=0.12, color="#2980B9")
    for _, row in df.iterrows():
        axes[0].annotate(
            f'{row["Silhouette score"]:.3f}',
            (row["k"], row["Silhouette score"]),
            textcoords="offset points", xytext=(0, 8),
            ha="center", fontsize=8, color="#2980B9",
            clip_on=False,
        )
    axes[0].set_xlabel("Number of clusters (k)", fontsize=11)
    axes[0].set_ylabel("Silhouette score", fontsize=11)
    axes[0].set_title("Internal validity\n(higher = better-defined clusters)", fontsize=10)
    axes[0].set_xticks([2, 3, 4, 5])

    if df["ARI vs diagnosis"].notna().any():
        axes[1].plot(df["k"], df["ARI vs diagnosis"], "o-",
                     color="#C1392B", lw=2, ms=8, zorder=3)
        axes[1].axhline(0, color="#AAAAAA", lw=0.8, ls="--")
        ari_min = df["ARI vs diagnosis"].min()
        ari_max = df["ARI vs diagnosis"].max()
        for _, row in df.iterrows():
            if row["ARI vs diagnosis"] is not None:
                a = row["ARI vs diagnosis"]
                if a == ari_min:
                    x_offset, y_offset, ha = 0, -14, "center"   # below dot
                elif a == ari_max:
                    x_offset, y_offset, ha = 12, 0, "left"      # right of dot
                else:
                    x_offset, y_offset, ha = 0, 10, "center"    # above dot
                axes[1].annotate(
                    f'{a:.3f}',
                    (row["k"], a),
                    textcoords="offset points", xytext=(x_offset, y_offset),
                    ha=ha, fontsize=8, color="#C1392B",
                    clip_on=False,
                )
        ymin, ymax = axes[1].get_ylim()
        axes[1].set_ylim(ymin - abs(ymin) * 0.8, ymax)
        axes[1].set_xlabel("Number of clusters (k)", fontsize=11)
        axes[1].set_ylabel("Adjusted Rand Index", fontsize=11)
        axes[1].set_title("External validity vs diagnosis\n(1.0 = perfect, 0 = random)",
                          fontsize=10)
        axes[1].set_xticks([2, 3, 4, 5])

    plt.tight_layout()
    fig.savefig(FIG_DIR / "fig5_cluster_validation.png", bbox_inches="tight")
    plt.close()

def main():
    print("=" * 55)
    print(" Step 3: Cluster Analysis")
    print("=" * 55)
    expr, meta, genes = load()
    plot_dendrogram(expr, meta)
    plot_pca(expr, meta)
    plot_gabaergic_coexpression(expr, genes, meta)
    plot_gene_heatmap(expr, meta, genes)
    plot_cluster_validation(expr, meta)
    print("\nAll figures saved to figures/")
    print("Run 04_pathway_analysis.py next.")


if __name__ == "__main__":
    main()