"""
04_pathway_analysis.py  (updated — all figures legible, no overlaps)
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
    print("\n── Differential expression ──────────────────────────────")
    scz_idx  = meta.index[meta["diagnosis"] == "Schizophrenia"]
    ctrl_idx = meta.index[meta["diagnosis"] == "Control"]

    rows = []
    for gene in raw.columns:
        scz_v  = raw.loc[raw.index.intersection(scz_idx),  gene].dropna()
        ctrl_v = raw.loc[raw.index.intersection(ctrl_idx), gene].dropna()
        if len(scz_v) < 3 or len(ctrl_v) < 3: continue
        t, p = stats.ttest_ind(scz_v, ctrl_v, equal_var=False)
        rows.append({"gene": gene,
                     "mean_scz":  round(scz_v.mean(), 4),
                     "mean_ctrl": round(ctrl_v.mean(), 4),
                     "log2fc":    round(scz_v.mean() - ctrl_v.mean(), 4),
                     "pval": p, "n_scz": len(scz_v), "n_ctrl": len(ctrl_v)})

    de = pd.DataFrame(rows).sort_values("pval").reset_index(drop=True)
    m  = len(de)
    de["rank"] = range(1, m + 1)
    de["padj"] = (de["pval"] * m / de["rank"]).clip(upper=1.0).round(4)
    de["significant"] = (de["padj"] < ALPHA) & (de["log2fc"].abs() >= MIN_FC)

    sig = de["significant"].sum()
    print(f"  Tested: {m}  |  Significant (FDR<{ALPHA}, |FC|≥{MIN_FC}): {sig}")
    if sig:
        print(f"  Hits: {de[de['significant']]['gene'].tolist()}")
    de.to_csv(OUT_DIR / "differential_expression.csv", index=False)
    return de


# ── Fig 6: Volcano ────────────────────────────────────────────────────────────
def plot_volcano(de, genes):
    print("\n── Fig 6: Volcano plot ──────────────────────────────────")
    g2s = dict(zip(genes["gene"], genes["subtype"]))

    fig, ax = plt.subplots(figsize=(9, 7))

    # Non-significant: small grey
    ns = de[~de["significant"]]
    ax.scatter(ns["log2fc"], -np.log10(ns["pval"] + 1e-10),
               c="#CCCCCC", s=25, alpha=0.55, zorder=2)

    # Significant: coloured, larger
    sig = de[de["significant"]]
    for _, row in sig.iterrows():
        color = SUBTYPE_COLORS.get(g2s.get(row["gene"], "SCZ candidate"), "#999")
        ax.scatter(row["log2fc"], -np.log10(row["pval"] + 1e-10),
                   c=color, s=90, alpha=0.95, zorder=4,
                   edgecolors="white", linewidths=0.7)

    # Manual label offsets tuned to avoid all overlaps
    LABEL_OFFSETS = {
        "NPY":    (-45, -18),
        "SST":    (-42,   6),
        "PVALB":  ( 10,   6),
        "PENK":   (-42,   6),
        "GAD1":   ( 10,  -14),
        "CALB2":  ( 10,   6),
        "CNR1":   ( 10,  -28),
        "PCP4":   (-42,  -14),
        "RASGRF2":(-48,   6),
        "SYNPR":  ( 10,   6),
        "TAC1":   (-42,  -14),
        "FEZ1":   ( 10,   6),
        "SNCG":   ( 10, -14),
        "NOS1AP": ( 10,   6),
        "CLDN5":  (-42,   6),
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

    # Reference lines
    y_thresh = -np.log10(ALPHA)
    ax.axhline(y_thresh, color="#888888", lw=0.9, ls="--", alpha=0.7)
    ax.axvline( MIN_FC, color="#CCCCCC", lw=0.8, ls=":", alpha=0.7)
    ax.axvline(-MIN_FC, color="#CCCCCC", lw=0.8, ls=":", alpha=0.7)
    ax.axvline(0, color="#AAAAAA", lw=0.8, alpha=0.4)

    # FDR label — placed inside plot area, not clipped
    x_range = de["log2fc"].max() - de["log2fc"].min()
    ax.text(de["log2fc"].min() + 0.03 * x_range, y_thresh + 0.12,
            f"FDR = {ALPHA}", fontsize=8, color="#888888", va="bottom")

    # Direction labels: bottom corners avoid collision with gene labels at top-left
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

    # Legend for subtypes present in significant genes
    present = {g2s.get(r["gene"], "") for _, r in sig.iterrows()}
    patches = [mpatches.Patch(color="#CCCCCC", label="Not significant")]
    patches += [mpatches.Patch(color=v, label=k)
                for k, v in SUBTYPE_COLORS.items() if k in present]
    ax.legend(handles=patches, title="Subtype", title_fontsize=8,
              frameon=True, edgecolor="#CCCCCC", fontsize=8,
              loc="lower left", bbox_to_anchor=(0.01, 0.12))

    fig.savefig(FIG_DIR / "fig6_volcano.png")
    plt.close()
    print(f"  {len(sig)} significant genes")


# ── Fig 7: Fold-change bar chart ──────────────────────────────────────────────
def plot_fc_barplot(de, genes):
    print("\n── Fig 7: FC bar chart ──────────────────────────────────")
    g2s = dict(zip(genes["gene"], genes["subtype"]))

    top = (de.reindex(de["log2fc"].abs().sort_values(ascending=False).index)
             .head(20)
             .sort_values("log2fc"))  # most negative at bottom

    colors     = [SCZ_COLOR if fc > 0 else CTRL_COLOR for fc in top["log2fc"]]
    edge_colors = ["black" if s else "none" for s in top["significant"]]

    # Wider right margin to fit p-values without overlap
    fig, ax = plt.subplots(figsize=(10, 9))
    fig.subplots_adjust(right=0.72)   # leave 28% right margin for p-values

    bars = ax.barh(
        range(len(top)), top["log2fc"].values,
        color=colors, alpha=0.82,
        edgecolor=edge_colors, linewidth=0.9,
        height=0.65
    )
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["gene"].values, fontsize=9)

    # Colour gene labels by subtype
    for tick, gene in zip(ax.get_yticklabels(), top["gene"].values):
        tick.set_color(SUBTYPE_COLORS.get(g2s.get(gene, ""), "#333333"))
        tick.set_fontweight("bold")

    ax.axvline(0, color="black", lw=0.9)
    ax.axvline( MIN_FC, color="#CCCCCC", lw=0.7, ls="--", alpha=0.7)
    ax.axvline(-MIN_FC, color="#CCCCCC", lw=0.7, ls="--", alpha=0.7)

    # Significance stars — to the right of the bar
    for i, (_, row) in enumerate(top.iterrows()):
        if row["significant"]:
            x_star = row["log2fc"] + (0.005 if row["log2fc"] > 0 else -0.005)
            ha_star = "left" if row["log2fc"] > 0 else "right"
            ax.text(x_star, i, " ★", va="center", ha=ha_star,
                    fontsize=11, color="black")

    # P-values in the right margin — well outside the bars
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

    # Direction legend
    dir_patches = [
        mpatches.Patch(color=SCZ_COLOR,  label="Higher in Schizophrenia"),
        mpatches.Patch(color=CTRL_COLOR, label="Lower in Schizophrenia"),
    ]
    # Subtype legend
    present = {g2s.get(g, "") for g in top["gene"].values}
    sub_patches = [mpatches.Patch(color=v, label=k)
                   for k, v in SUBTYPE_COLORS.items()
                   if k in present and k != "SCZ candidate"]

    # Two separate legends — direction (inside) and subtype (outside right)
    ax.legend(handles=dir_patches, loc="lower right",
              frameon=True, edgecolor="#CCCCCC", fontsize=8)
    ax.figure.legend(
        handles=sub_patches,
        loc="upper right", bbox_to_anchor=(1.0, 0.88),
        title="Subtype", title_fontsize=8,
        frameon=True, edgecolor="#CCCCCC", fontsize=8
    )

    fig.savefig(FIG_DIR / "fig7_barplot_fc.png", bbox_inches="tight")
    plt.close()
    print("  Saved → fig7_barplot_fc.png")


# ── Fig 8: Pathway enrichment ─────────────────────────────────────────────────
def plot_pathway_enrichment(de, genes):
    print("\n── Fig 8: Pathway enrichment ────────────────────────────")
    g2s   = dict(zip(genes["gene"], genes["subtype"]))
    all_g = set(de["gene"])
    deg_g = set(de[de["significant"]]["gene"])

    rows = []
    for sub in genes["subtype"].dropna().unique():
        in_set = {g for g, s in g2s.items() if s == sub} & all_g
        if len(in_set) < 2: continue
        a = len(deg_g & in_set)
        b = len(deg_g - in_set)
        c = len((all_g - deg_g) & in_set)
        d = len((all_g - deg_g) - in_set)
        _, p = fisher_exact([[a, b],[c, d]], alternative="greater")
        rows.append({"subtype": sub, "pct_deg": a/len(in_set)*100,
                     "n_deg": a, "n_genes": len(in_set), "pval": p})

    df = pd.DataFrame(rows).sort_values("pct_deg", ascending=True)

    # Taller figure so text doesn't overlap
    fig, ax = plt.subplots(figsize=(9, max(5, len(df) * 0.7)))
    fig.subplots_adjust(right=0.70)   # right margin for labels

    colors = [SUBTYPE_COLORS.get(s, "#CCCCCC") for s in df["subtype"]]
    y_pos  = range(len(df))
    bars   = ax.barh(y_pos, df["pct_deg"].values,
                     color=colors, alpha=0.85, height=0.6)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df["subtype"].values, fontsize=11)
    # Colour subtype labels to match bars
    for tick, sub in zip(ax.get_yticklabels(), df["subtype"].values):
        tick.set_color(SUBTYPE_COLORS.get(sub, "#333"))
        tick.set_fontweight("bold")

    # Fraction and p-value labels — outside bars in right margin
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

    # Reference line at the background DEG rate
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
    print(" Step 4: Pathway Analysis (updated figures)")
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
