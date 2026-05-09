"""
04_pathway_analysis.py
Differential expression and pathway enrichment analysis.

Reads preprocessed log2 expression data from 02_preprocess.py and
produces three figures and two output tables:

  fig6_volcano.png          Volcano plot of all genes. Significant genes
                            colored by interneuron subtype.

  fig7_barplot_fc.png       Horizontal bar chart of top 20 genes by |log2FC|,
                            with significance stars and p-values.

  fig8_pathway_enrichment.png
                            Over-representation analysis: % of each interneuron
                            subtype's genes that are DEGs. Fisher's exact test.

  output/differential_expression.csv   Full DE results for all genes.
  output/top_genes.csv                 Top 10 genes by |log2FC|.

Statistical methods:
Differential expression: Welch's t-test (unequal variance assumed, appropriate
for unequal group sizes) per gene, followed by Benjamini-Hochberg FDR correction.

Thresholds: FDR < 0.10, |log2FC| >= 0.10. These thresholds are intentionally
relaxed relative to genome-wide standards (FDR < 0.05, |FC| > 1.5) because:
  1. This is a curated 48-gene panel, not an unbiased genome-wide screen.
     The multiple testing burden is far lower (48 tests, not ~20,000).
  2. Effect sizes for interneuron markers in postmortem schizophrenia are
     consistently modest (log2FC 0.3-0.8 in the literature).
  3. Sample size (n=34 DLPFC donors) limits power to detect small effects at strict thresholds.
All thresholds are disclosed and justified.

Pathway enrichment: Over-representation analysis (ORA) using Fisher's exact
test. Gene sets are the manually curated interneuron subtype labels defined
in 02_preprocess.py. Enrichment against external databases (KEGG, GO, Reactome)
was not implemented in this pipeline.

References:
Guillozet-Bongaarts AL, et al. (2014). Altered gene expression in the
  dorsolateral prefrontal cortex of individuals with schizophrenia.
  Molecular Psychiatry, 19(4):478-485. doi:10.1038/mp.2013.30. PMID: 23528911.

Khatri P, Sirota M, Butte AJ (2012). Ten years of pathway analysis: current
  approaches and outstanding challenges. PLoS Computational Biology,
  8(2):e1002375. doi:10.1371/journal.pcbi.1002375.

Benjamini Y, Hochberg Y (1995). Controlling the false discovery rate: a
  practical and powerful approach to multiple testing. J R Stat Soc B,
  57:289-300.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.stats import fisher_exact
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

DATA_DIR = Path(__file__).parent.parent / "data"
FIG_DIR  = Path(__file__).parent.parent / "figures"
OUT_DIR  = Path(__file__).parent.parent / "output"
FIG_DIR.mkdir(exist_ok=True)
OUT_DIR.mkdir(exist_ok=True)

ALPHA   = 0.10
MIN_FC  = 0.10

SCZ_COLOR  = "#C1392B"
CTRL_COLOR = "#2980B9"

SUBTYPE_COLORS = {
    "PV+":             "#7B2D8B",
    "SST+":            "#2471A3",
    "VIP+":            "#1E8449",
    "CB+":             "#D4AC0D",
    "Pan-GABA":        "#BA4A00",
    "Excitatory":      "#7D3C98",
    "Oligodendrocyte": "#566573",
    "Risk/Signaling":  "#909497",
    "Synaptic":        "#B2BABB",
    "SCZ candidate":   "#D5D8DC",
}

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor":   "white",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.linewidth":    0.8,
    "font.family":       "DejaVu Sans",
    "savefig.dpi":       180,
    "savefig.bbox":      "tight",
    "savefig.facecolor": "white",
})


def load():
    raw   = pd.read_csv(DATA_DIR / "expression_matrix_raw.csv", index_col=0)
    meta  = pd.read_csv(DATA_DIR / "donor_metadata.csv",        index_col=0)
    genes = pd.read_csv(DATA_DIR / "gene_annotations.csv")
    return raw, meta.reindex(raw.index), genes


def differential_expression(raw, meta):
    """
    Welch's t-test per gene + Benjamini-Hochberg FDR correction.
    Returns DataFrame with log2FC, p-value, adjusted p-value, and significance flag.
    """
    print("\nDifferential expression")

    scz_idx  = meta.index[meta["diagnosis"] == "Schizophrenia"]
    ctrl_idx = meta.index[meta["diagnosis"] == "Control"]

    print("  SCZ: " + str(len(scz_idx)) + "  |  Control: " + str(len(ctrl_idx)))

    rows = []
    for gene in raw.columns:
        scz_v  = raw.loc[raw.index.intersection(scz_idx),  gene].dropna()
        ctrl_v = raw.loc[raw.index.intersection(ctrl_idx), gene].dropna()

        # Need at least 3 donors per group to run a t-test.
        if len(scz_v) >= 3 and len(ctrl_v) >= 3:
            t, p = stats.ttest_ind(scz_v, ctrl_v, equal_var=False)
            rows.append({
                "gene":      gene,
                "mean_scz":  round(scz_v.mean(), 4),
                "mean_ctrl": round(ctrl_v.mean(), 4),
                "log2fc":    round(scz_v.mean() - ctrl_v.mean(), 4),
                "pval":      p,
                "n_scz":     len(scz_v),
                "n_ctrl":    len(ctrl_v),
            })

    de = pd.DataFrame(rows).sort_values("pval").reset_index(drop=True)
    m  = len(de)

    de["rank"] = range(1, m + 1)

    # Benjamini-Hochberg FDR: multiply each sorted p-value by m/rank,
    # clip at 1, then enforce monotonicity from the largest rank downward
    # so that padj[i] <= padj[i+1] (cummin on the reversed series).
    de["padj"]        = (de["pval"] * m / de["rank"]).clip(upper=1.0)
    de["padj"]        = de["padj"][::-1].cummin()[::-1].round(4)
    de["significant"] = (de["padj"] < ALPHA) & (de["log2fc"].abs() >= MIN_FC)

    sig = de["significant"].sum()
    print("  Genes tested: " + str(m) +
          "  |  Significant (FDR<" + str(ALPHA) +
          ", |FC|>=" + str(MIN_FC) + "): " + str(sig))

    if sig:
        print("  Hits: " + str(de[de["significant"]]["gene"].tolist()))

    de.to_csv(OUT_DIR / "differential_expression.csv", index=False)
    return de


def plot_volcano(de, genes):
    """
    Volcano plot: log2FC (x) vs -log10(p-value) (y).
    Significant genes labeled and colored by interneuron subtype.
    """
    print("\nFig 6: Volcano plot")

    g2s = dict(zip(genes["gene"], genes["subtype"]))

    fig, ax = plt.subplots(figsize=(9, 7))

    # Non-significant genes: small gray dots.
    ns = de[~de["significant"]]
    ax.scatter(ns["log2fc"], -np.log10(ns["pval"] + 1e-10),
               c="#CCCCCC", s=25, alpha=0.55, zorder=2)

    # Significant genes: larger dots colored by subtype.
    sig = de[de["significant"]]
    for i in range(len(sig)):
        row   = sig.iloc[i]
        color = SUBTYPE_COLORS.get(g2s.get(row["gene"], "SCZ candidate"), "#999")
        ax.scatter(row["log2fc"], -np.log10(row["pval"] + 1e-10),
                   c=color, s=90, alpha=0.95, zorder=4,
                   edgecolors="white", linewidths=0.7)

    # Manual label offsets tuned to avoid overlaps.
    LABEL_OFFSETS = {
        "NPY":     (-45, -18),
        "SST":     (-42,   6),
        "PVALB":   ( 10,   6),
        "PENK":    (-42,   6),
        "GAD1":    ( 10, -14),
        "CALB2":   ( 22,   6),
        "CNR1":    (-42, -28),
        "PCP4":    (-42,   6),
        "RASGRF2": (-55, -14),
        "SYNPR":   ( 10,   6),
        "TAC1":    (-42, -14),
        "FEZ1":    ( 10,   6),
        "SNCG":    ( 10, -14),
        "NOS1AP":  ( 10,   6),
        "CLDN5":   (-42,   6),
    }

    for i in range(len(sig)):
        row   = sig.iloc[i]
        x0    = row["log2fc"]
        y0    = -np.log10(row["pval"] + 1e-10)
        color = SUBTYPE_COLORS.get(g2s.get(row["gene"], "SCZ candidate"), "#333")
        xoff, yoff = LABEL_OFFSETS.get(row["gene"], (10, 6))
        ax.annotate(
            row["gene"],
            xy=(x0, y0),
            xytext=(xoff, yoff), textcoords="offset points",
            fontsize=9, fontweight="bold", color=color,
            arrowprops=dict(arrowstyle="-", color="#BBBBBB", lw=0.6),
        )

    y_thresh = -np.log10(ALPHA)
    ax.axhline(y_thresh, color="#888888", lw=0.9, ls="--", alpha=0.7)
    ax.axvline( MIN_FC, color="#CCCCCC", lw=0.8, ls=":", alpha=0.7)
    ax.axvline(-MIN_FC, color="#CCCCCC", lw=0.8, ls=":", alpha=0.7)
    ax.axvline(0,       color="#AAAAAA", lw=0.8, alpha=0.4)

    x_range = de["log2fc"].max() - de["log2fc"].min()
    ax.text(de["log2fc"].min() + 0.03 * x_range, y_thresh + 0.12,
            "FDR = " + str(ALPHA), fontsize=8, color="#888888", va="bottom")

    ax.text(0.98, 0.06, "Higher in SCZ", transform=ax.transAxes,
            fontsize=9, color=SCZ_COLOR, ha="right", va="bottom")
    ax.text(0.02, 0.06, "Lower in SCZ", transform=ax.transAxes,
            fontsize=9, color=CTRL_COLOR, ha="left", va="bottom")

    ax.set_xlabel("log2 Fold Change  (Schizophrenia - Control)", fontsize=11)
    ax.set_ylabel("-log10(p-value)", fontsize=11)
    ax.set_title(
        "Differential Expression - DLPFC, Schizophrenia vs Control\n"
        "FDR < " + str(ALPHA) +
        "  |  |log2FC| >= " + str(MIN_FC) +
        "  |  Color = interneuron subtype",
        fontsize=11
    )

    # Build the set of subtypes present among significant genes.
    present = set()
    for i in range(len(sig)):
        present.add(g2s.get(sig.iloc[i]["gene"], ""))

    patches = [mpatches.Patch(color="#CCCCCC", label="Not significant")]
    for k in SUBTYPE_COLORS:
        if k in present:
            patches.append(mpatches.Patch(color=SUBTYPE_COLORS[k], label=k))

    ax.legend(handles=patches, title="Subtype", title_fontsize=8,
              frameon=True, edgecolor="#CCCCCC", fontsize=8,
              loc="upper right")

    fig.savefig(FIG_DIR / "fig6_volcano.png")
    plt.close()
    print("  " + str(len(sig)) + " significant genes labelled")


def plot_fc_barplot(de, genes):
    """
    Top 20 genes by |log2FC|, sorted ascending. Significant genes marked with star.
    P-values printed in right margin. Color = direction (red = up in SCZ).
    """
    print("\nFig 7: FC bar chart")

    g2s = dict(zip(genes["gene"], genes["subtype"]))

    top = (de.reindex(de["log2fc"].abs().sort_values(ascending=False).index)
             .head(20)
             .sort_values("log2fc"))

    # Build bar colors and edge colors as explicit lists.
    colors      = []
    edge_colors = []
    for fc in top["log2fc"]:
        if fc > 0:
            colors.append(SCZ_COLOR)
        else:
            colors.append(CTRL_COLOR)
    for s in top["significant"]:
        if s:
            edge_colors.append("black")
        else:
            edge_colors.append("none")

    fig, ax = plt.subplots(figsize=(10, 9))
    fig.subplots_adjust(right=0.72)

    ax.barh(
        range(len(top)), top["log2fc"].values,
        color=colors, alpha=0.82,
        edgecolor=edge_colors, linewidth=0.9,
        height=0.65
    )

    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["gene"].values, fontsize=9)

    for tick in ax.get_yticklabels():
        gene = tick.get_text()
        tick.set_color(SUBTYPE_COLORS.get(g2s.get(gene, ""), "#333333"))
        tick.set_fontweight("bold")

    ax.axvline(0,       color="black",   lw=0.9)
    ax.axvline( MIN_FC, color="#CCCCCC", lw=0.7, ls="--", alpha=0.7)
    ax.axvline(-MIN_FC, color="#CCCCCC", lw=0.7, ls="--", alpha=0.7)

    # Significance stars next to bars.
    for i in range(len(top)):
        row = top.iloc[i]
        if row["significant"]:
            if row["log2fc"] > 0:
                x_star = row["log2fc"] + 0.005
                ha_star = "left"
            else:
                x_star = row["log2fc"] - 0.005
                ha_star = "right"
            ax.text(x_star, i, " *", va="center", ha=ha_star,
                    fontsize=11, color="black")

    # P-value annotations in the right margin.
    x_pval = ax.get_xlim()[1] * 1.05
    for i in range(len(top)):
        row = top.iloc[i]
        if row["significant"]:
            style = "bold"
            color = "#333333"
        else:
            style = "normal"
            color = "#AAAAAA"
        ax.text(x_pval, i, "p = " + "%.3f" % row["pval"],
                va="center", ha="left", fontsize=8,
                fontweight=style, color=color,
                transform=ax.transData, clip_on=False)

    ax.set_xlabel("log2 Fold Change  (Schizophrenia - Control)", fontsize=11)
    ax.set_title(
        "Top 20 Genes by Effect Size - DLPFC\n"
        "* = significant after FDR correction  |  Gene color = subtype",
        fontsize=11
    )

    # Build the set of subtypes present in the top genes.
    present = set()
    for g in top["gene"].values:
        present.add(g2s.get(g, ""))

    sub_patches = []
    for k in SUBTYPE_COLORS:
        if k in present and k != "SCZ candidate":
            sub_patches.append(mpatches.Patch(color=SUBTYPE_COLORS[k], label=k))

    dir_patches = [
        mpatches.Patch(color=SCZ_COLOR,  label="Higher in Schizophrenia"),
        mpatches.Patch(color=CTRL_COLOR, label="Lower in Schizophrenia"),
    ]

    ax.figure.legend(
        handles=sub_patches,
        loc="upper right", bbox_to_anchor=(1.0, 0.88),
        title="Subtype", title_fontsize=8,
        frameon=True, edgecolor="#CCCCCC", fontsize=8
    )
    ax.figure.legend(
        handles=dir_patches,
        loc="upper right", bbox_to_anchor=(1.0, 0.70),
        title="Direction", title_fontsize=8,
        frameon=True, edgecolor="#CCCCCC", fontsize=8
    )

    fig.savefig(FIG_DIR / "fig7_barplot_fc.png", bbox_inches="tight")
    plt.close()
    print("  Saved -> fig7_barplot_fc.png")


def plot_pathway_enrichment(de, genes):
    """
    Over-representation analysis (Fisher's exact test) per interneuron subtype.
    Tests whether differentially expressed genes are enriched within each
    manually curated subtype category.
    Background DEG rate shown as a dashed reference line.
    """
    print("\nFig 8: Pathway enrichment")

    g2s   = dict(zip(genes["gene"], genes["subtype"]))
    all_g = set(de["gene"])
    deg_g = set(de[de["significant"]]["gene"])

    rows = []
    for sub in genes["subtype"].dropna().unique():

        # Build the set of genes in this subtype that are in our tested panel.
        in_set = set()
        for g in g2s:
            if g2s[g] == sub and g in all_g:
                in_set.add(g)

        if len(in_set) < 2:
            pass  # skip subtypes with fewer than 2 genes in the panel
        else:
            # 2x2 contingency table for Fisher's exact test:
            #   a = DEG and in subtype
            #   b = DEG and not in subtype
            #   c = not DEG and in subtype
            #   d = not DEG and not in subtype
            a = len(deg_g & in_set)
            b = len(deg_g - in_set)
            c = len((all_g - deg_g) & in_set)
            d = len((all_g - deg_g) - in_set)

            _, p = fisher_exact([[a, b], [c, d]], alternative="greater")

            rows.append({
                "subtype":  sub,
                "pct_deg":  a / len(in_set) * 100,
                "n_deg":    a,
                "n_genes":  len(in_set),
                "pval":     p,
            })

    df = pd.DataFrame(rows).sort_values("pct_deg", ascending=True)

    fig, ax = plt.subplots(figsize=(9, max(5, len(df) * 0.7)))
    fig.subplots_adjust(right=0.70)

    # Build bar color list for each subtype row.
    colors = []
    for s in df["subtype"]:
        colors.append(SUBTYPE_COLORS.get(s, "#CCCCCC"))

    y_pos = range(len(df))
    ax.barh(y_pos, df["pct_deg"].values, color=colors, alpha=0.85, height=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["subtype"].values, fontsize=11)

    for tick in ax.get_yticklabels():
        sub = tick.get_text()
        tick.set_color(SUBTYPE_COLORS.get(sub, "#333"))
        tick.set_fontweight("bold")

    x_max = max(df["pct_deg"].max(), 10)

    # Annotate each bar with gene counts and p-value.
    for i in range(len(df)):
        row = df.iloc[i]

        if row["pval"] < 0.05:
            sig_star = " *"
            color    = "#333333"
            weight   = "bold"
        else:
            sig_star = ""
            color    = "#888888"
            weight   = "normal"

        label = (str(row["n_deg"]) + "/" + str(row["n_genes"]) +
                 " genes   p = " + "%.3f" % row["pval"] + sig_star)

        ax.text(
            x_max * 1.04, i, label,
            va="center", ha="left", fontsize=9,
            color=color, fontweight=weight,
            clip_on=False
        )

    ax.set_xlim(0, x_max)
    ax.set_xlabel("% of subtype genes that are differentially expressed", fontsize=11)
    ax.set_title(
        "Pathway Enrichment by Interneuron Subtype\n"
        "Over-Representation Analysis (Fisher's exact test)  |  * p < 0.05",
        fontsize=11
    )

    if len(all_g) > 0:
        bg_rate = len(deg_g) / len(all_g) * 100
    else:
        bg_rate = 0

    ax.axvline(bg_rate, color="#888888", lw=0.9, ls="--", alpha=0.7)
    ax.text(bg_rate + 0.5, len(df) - 0.3,
            "Background\n(" + "%.0f" % bg_rate + "%)",
            fontsize=7.5, color="#888888", va="top")

    fig.savefig(FIG_DIR / "fig8_pathway_enrichment.png", bbox_inches="tight")
    plt.close()
    print("  Saved -> fig8_pathway_enrichment.png")


def print_summary(de):
    print("\nTop genes by |log2FC|")
    top = (de.reindex(de["log2fc"].abs().sort_values(ascending=False).index)
             .head(10))
    print(top[["gene", "log2fc", "pval", "padj", "significant"]].to_string(index=False))
    top.to_csv(OUT_DIR / "top_genes.csv", index=False)


def main():
    print("=" * 55)
    print(" Step 4: Pathway Analysis")
    print("=" * 55)

    raw, meta, genes = load()
    de = differential_expression(raw, meta)

    plot_volcano(de, genes)
    plot_fc_barplot(de, genes)
    plot_pathway_enrichment(de, genes)
    print_summary(de)

    print("\nAnalysis complete.")


if __name__ == "__main__":
    main()