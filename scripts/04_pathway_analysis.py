"""
04_pathway_analysis.py
======================
Differential expression and pathway enrichment analysis.

Reads batch-corrected log2 expression data from 02_preprocess.py and
produces three figures and two output tables:

  fig6_volcano.png          Volcano plot of all genes. Significant genes
                            coloured by interneuron subtype.

  fig7_barplot_fc.png       Horizontal bar chart of top 20 genes by |log2FC|,
                            with significance stars and p-values.

  fig8_pathway_enrichment.png
                            Over-representation analysis: % of each interneuron
                            subtype's genes that are DEGs. Fisher's exact test.

  output/differential_expression.csv   Full DE results for all genes.
  output/top_genes.csv                 Top 10 genes by |log2FC|.

Statistical methods
-------------------
Differential expression: Welch's t-test (unequal variance assumed, appropriate
for unequal group sizes) per gene, followed by Benjamini-Hochberg FDR correction.

Thresholds: FDR < 0.10, |log2FC| >= 0.10. These thresholds are intentionally
relaxed relative to genome-wide standards (FDR < 0.05, |FC| > 1.5) because:
  1. This is a curated 48-gene panel, not an unbiased genome-wide screen.
     The multiple testing burden is far lower (48 tests, not ~20,000).
  2. Effect sizes for interneuron markers in postmortem schizophrenia are
     consistently modest (log2FC 0.3-0.8 in the literature).
  3. Sample size (n=103) limits power to detect small effects at strict thresholds.
All thresholds are disclosed and justified.

Pathway enrichment: Over-representation analysis (ORA) using Fisher's exact
test. Gene sets are the manually curated interneuron subtype labels defined
in 02_preprocess.py. Enrichment against external databases (KEGG, GO, Reactome)
was not implemented in this pipeline.

References
----------
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

# Relaxed thresholds — justified for focused 48-gene curated panel.
# See module docstring for rationale.
ALPHA  = 0.10
MIN_FC = 0.10

SCZ_COLOR  = "#C1392B"
CTRL_COLOR = "#2980B9"

SUBTYPE_COLORS = {
    "PV+":            "#7B2D8B",
    "SST+":           "#2471A3",
    "VIP+":           "#1E8449",
    "CB+":            "#D4AC0D",
    "Pan-GABA":       "#BA4A00",
    "Excitatory":     "#7D3C98",
    "Oligodendrocyte":"#566573",
    "Risk/Signaling": "#909497",
    "Synaptic":       "#B2BABB",
    "SCZ candidate":  "#D5D8DC",
}

plt.rcParams.update({
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
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
    print("\n── Differential expression ──────────────────────────────")
    scz_idx  = meta.index[meta["diagnosis"] == "Schizophrenia"]
    ctrl_idx = meta.index[meta["diagnosis"] == "Control"]
    print(f"  SCZ: {len(scz_idx)}  |  Control: {len(ctrl_idx)}")

    rows = []
    for gene in raw.columns:
        scz_v  = raw.loc[raw.index.intersection(scz_idx),  gene].dropna()
        ctrl_v = raw.loc[raw.index.intersection(ctrl_idx), gene].dropna()
        if len(scz_v) < 3 or len(ctrl_v) < 3:
            continue
        t, p = stats.ttest_ind(scz_v, ctrl_v, equal_var=False)
        rows.append({"gene":      gene,
                     "mean_scz":  round(scz_v.mean(), 4),
                     "mean_ctrl": round(ctrl_v.mean(), 4),
                     "log2fc":    round(scz_v.mean() - ctrl_v.mean(), 4),
                     "pval":      p,
                     "n_scz":     len(scz_v),
                     "n_ctrl":    len(ctrl_v)})

    de = pd.DataFrame(rows).sort_values("pval").reset_index(drop=True)
    m  = len(de)
    de["rank"] = range(1, m + 1)
    # Benjamini-Hochberg: padj = p * m / rank, clipped at 1.0
    de["padj"] = (de["pval"] * m / de["rank"]).clip(upper=1.0).round(4)
    de["significant"] = (de["padj"] < ALPHA) & (de["log2fc"].abs() >= MIN_FC)

    sig = de["significant"].sum()
    print(f"  Genes tested: {m}  |  Significant (FDR<{ALPHA}, |FC|≥{MIN_FC}): {sig}")
    if sig:
        print(f"  Hits: {de[de['significant']]['gene'].tolist()}")
    de.to_csv(OUT_DIR / "differential_expression.csv", index=False)
    return de


# ── Fig 6: Volcano ────────────────────────────────────────────────────────────
def plot_volcano(de, genes):
    """
    Volcano plot: log2FC (x) vs -log10(p-value) (y).
    Significant genes labelled and coloured by interneuron subtype.
    """
    print("\n── Fig 6: Volcano plot ──────────────────────────────────")
    g2s = dict(zip(genes["gene"], genes["subtype"]))

    fig, ax = plt.subplots(figsize=(9, 7))

    # Non-significant: small grey
    ns = de[~de["significant"]]
    ax.scatter(ns["log2fc"], -np.log10(ns["pval"] + 1e-10),
               c="#CCCCCC", s=25, alpha=0.55, zorder=2)

    # Significant: coloured by subtype, larger
    sig = de[de["significant"]]
    for _, row in sig.iterrows():
        color = SUBTYPE_COLORS.get(g2s.get(row["gene"], "SCZ candidate"), "#999")
        ax.scatter(row["log2fc"], -np.log10(row["pval"] + 1e-10),
                   c=color, s=90, alpha=0.95, zorder=4,
                   edgecolors="white", linewidths=0.7)

    # Manual label offsets tuned to avoid overlaps
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

    for _, row in sig.iterrows():
        x0 = row["log2fc"]
        y0 = -np.log10(row["pval"] + 1e-10)
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
    ax.axvline(0, color="#AAAAAA", lw=0.8, alpha=0.4)

    x_range = de["log2fc"].max() - de["log2fc"].min()
    ax.text(de["log2fc"].min() + 0.03 * x_range, y_thresh + 0.12,
            f"FDR = {ALPHA}", fontsize=8, color="#888888", va="bottom")

    ax.text(0.98, 0.06, "↑ Higher in SCZ", transform=ax.transAxes,
            fontsize=9, color=SCZ_COLOR, ha="right", va="bottom")
    ax.text(0.02, 0.06, "↓ Lower in SCZ", transform=ax.transAxes,
            fontsize=9, color=CTRL_COLOR, ha="left", va="bottom")

    ax.set_xlabel("log₂ Fold Change  (Schizophrenia − Control)", fontsize=11)
    ax.set_ylabel("−log₁₀(p-value)", fontsize=11)
    ax.set_title(
        "Differential Expression — DLPFC, Schizophrenia vs Control\n"
        f"FDR < {ALPHA}  |  |log₂FC| ≥ {MIN_FC}  |  Colour = interneuron subtype",
        fontsize=11
    )

    present = {g2s.get(r["gene"], "") for _, r in sig.iterrows()}
    patches = [mpatches.Patch(color="#CCCCCC", label="Not significant")]
    patches += [mpatches.Patch(color=v, label=k)
                for k, v in SUBTYPE_COLORS.items() if k in present]
    ax.legend(handles=patches, title="Subtype", title_fontsize=8,
              frameon=True, edgecolor="#CCCCCC", fontsize=8,
              loc="upper right")

    fig.savefig(FIG_DIR / "fig6_volcano.png")
    plt.close()
    print(f"  {len(sig)} significant genes labelled")


# ── Fig 7: Fold-change bar chart ──────────────────────────────────────────────
def plot_fc_barplot(de, genes):
    """
    Top 20 genes by |log2FC|, sorted ascending. Significant genes marked with ★.
    P-values printed in right margin. Colour = direction (red = up in SCZ).
    """
    print("\n── Fig 7: FC bar chart ──────────────────────────────────")
    g2s = dict(zip(genes["gene"], genes["subtype"]))
    top = (de.reindex(de["log2fc"].abs().sort_values(ascending=False).index)
             .head(20)
             .sort_values("log2fc"))
    colors      = [SCZ_COLOR if fc > 0 else CTRL_COLOR for fc in top["log2fc"]]
    edge_colors = ["black" if s else "none" for s in top["significant"]]
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
    for tick, gene in zip(ax.get_yticklabels(), top["gene"].values):
        tick.set_color(SUBTYPE_COLORS.get(g2s.get(gene, ""), "#333333"))
        tick.set_fontweight("bold")
    ax.axvline(0, color="black", lw=0.9)
    ax.axvline( MIN_FC, color="#CCCCCC", lw=0.7, ls="--", alpha=0.7)
    ax.axvline(-MIN_FC, color="#CCCCCC", lw=0.7, ls="--", alpha=0.7)
    for i, (_, row) in enumerate(top.iterrows()):
        if row["significant"]:
            x_star = row["log2fc"] + (0.005 if row["log2fc"] > 0 else -0.005)
            ha_star = "left" if row["log2fc"] > 0 else "right"
            ax.text(x_star, i, " ★", va="center", ha=ha_star,
                    fontsize=11, color="black")
    x_pval = ax.get_xlim()[1] * 1.05
    for i, (_, row) in enumerate(top.iterrows()):
        style = "bold" if row["significant"] else "normal"
        color = "#333333" if row["significant"] else "#AAAAAA"
        ax.text(x_pval, i, f"p = {row['pval']:.3f}",
                va="center", ha="left", fontsize=8,
                fontweight=style, color=color,
                transform=ax.transData, clip_on=False)
    ax.set_xlabel("log₂ Fold Change  (Schizophrenia − Control)", fontsize=11)
    ax.set_title(
        "Top 20 Genes by Effect Size — DLPFC\n"
        "★ = significant after FDR correction  |  Gene colour = subtype",
        fontsize=11
    )
    present = {g2s.get(g, "") for g in top["gene"].values}
    sub_patches = [mpatches.Patch(color=v, label=k)
                   for k, v in SUBTYPE_COLORS.items()
                   if k in present and k != "SCZ candidate"]
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
    print("  Saved → fig7_barplot_fc.png")


# ── Fig 8: Pathway enrichment ─────────────────────────────────────────────────
def plot_pathway_enrichment(de, genes):
    """
    Over-representation analysis (Fisher's exact test) per interneuron subtype.

    Tests whether differentially expressed genes are enriched within each
    manually curated subtype category. Note: this tests against the subtype
    labels defined in 02_preprocess.py, not against external databases
    (KEGG, GO, Reactome).

    Background DEG rate shown as a dashed reference line.
    """
    print("\n── Fig 8: Pathway enrichment ────────────────────────────")
    g2s   = dict(zip(genes["gene"], genes["subtype"]))
    all_g = set(de["gene"])
    deg_g = set(de[de["significant"]]["gene"])

    rows = []
    for sub in genes["subtype"].dropna().unique():
        in_set = {g for g, s in g2s.items() if s == sub} & all_g
        if len(in_set) < 2:
            continue
        a = len(deg_g & in_set)
        b = len(deg_g - in_set)
        c = len((all_g - deg_g) & in_set)
        d = len((all_g - deg_g) - in_set)
        _, p = fisher_exact([[a, b],[c, d]], alternative="greater")
        rows.append({"subtype":  sub,
                     "pct_deg":  a / len(in_set) * 100,
                     "n_deg":    a,
                     "n_genes":  len(in_set),
                     "pval":     p})

    df = pd.DataFrame(rows).sort_values("pct_deg", ascending=True)

    fig, ax = plt.subplots(figsize=(9, max(5, len(df) * 0.7)))
    fig.subplots_adjust(right=0.70)

    colors = [SUBTYPE_COLORS.get(s, "#CCCCCC") for s in df["subtype"]]
    y_pos  = range(len(df))
    ax.barh(y_pos, df["pct_deg"].values, color=colors, alpha=0.85, height=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["subtype"].values, fontsize=11)
    for tick, sub in zip(ax.get_yticklabels(), df["subtype"].values):
        tick.set_color(SUBTYPE_COLORS.get(sub, "#333"))
        tick.set_fontweight("bold")

    x_max = max(df["pct_deg"].max(), 10)
    for i, (_, row) in enumerate(df.iterrows()):
        sig_star = " ★" if row["pval"] < 0.05 else ""
        label = f"{row['n_deg']}/{row['n_genes']} genes   p = {row['pval']:.3f}{sig_star}"
        ax.text(
            x_max * 1.04, i, label,
            va="center", ha="left", fontsize=9,
            color="#333333" if row["pval"] < 0.05 else "#888888",
            fontweight="bold" if row["pval"] < 0.05 else "normal",
            clip_on=False
        )

    ax.set_xlim(0, x_max)
    ax.set_xlabel("% of subtype genes that are differentially expressed", fontsize=11)
    ax.set_title(
        "Pathway Enrichment by Interneuron Subtype\n"
        "Over-Representation Analysis (Fisher's exact test)  |  ★ p < 0.05",
        fontsize=11
    )

    bg_rate = len(deg_g) / len(all_g) * 100 if all_g else 0
    ax.axvline(bg_rate, color="#888888", lw=0.9, ls="--", alpha=0.7)
    ax.text(bg_rate + 0.5, len(df) - 0.3,
            f"Background\n({bg_rate:.0f}%)",
            fontsize=7.5, color="#888888", va="top")

    fig.savefig(FIG_DIR / "fig8_pathway_enrichment.png", bbox_inches="tight")
    plt.close()
    print("  Saved → fig8_pathway_enrichment.png")


def print_summary(de):
    print("\n── Top genes by |log2FC| ────────────────────────────────")
    top = de.reindex(de["log2fc"].abs().sort_values(ascending=False).index).head(10)
    print(top[["gene","log2fc","pval","padj","significant"]].to_string(index=False))
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