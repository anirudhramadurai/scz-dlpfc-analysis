"""
02_preprocess.py
Filters, normalizes, and annotates GEO expression data for DLPFC-only analysis.

Key steps:
1. Filter to Schizophrenia and Control DLPFC donors only
   (GSE53987 also contains bipolar disorder and MDD donors, and hippocampus
   and striatum samples, all excluded here)
2. Log2 transform expression values
3. Detect processing batch from GSM ID prefix; skip correction if single cohort
4. Aggregate to mean per donor x gene
5. Z-score normalize across donors per gene
6. Annotate genes with interneuron subtype labels

Batch correction note:
After filtering to DLPFC-only samples (n=34), all donors fall within the
GSM1304xxx processing cohort. The correct_batch() function detects this
automatically and skips correction when only one batch is present. The
42% PC1 variance previously attributed to a batch effect was observed in
the full 205-sample dataset and reflected biological differences between
brain regions (DLPFC, hippocampus, striatum), not a technical artifact.
No batch correction is applied to the DLPFC-only data. The correct_batch()
function is retained for use with multi-region data.

Note: explicit covariates (RIN, PMI, medication) were not available in this
dataset.

References:
Guillozet-Bongaarts AL, et al. (2014). Altered gene expression in the
  dorsolateral prefrontal cortex of individuals with schizophrenia.
  Molecular Psychiatry, 19(4):478-485. doi:10.1038/mp.2013.30. PMID: 23528911.

Lanz TA, et al. (2019). Postmortem transcriptional profiling reveals widespread
  increase in inflammation in schizophrenia. Translational Psychiatry, 9(1):151.
  doi:10.1038/s41398-019-0492-8. PMID: 31123247. GEO: GSE53987.

Dienel SJ, Lewis DA (2019). Alterations in cortical interneurons and cognitive
  function in schizophrenia. Neurobiology of Disease, 131:104208.
  PMC6309598.

Leek JT, et al. (2010). Tackling the widespread and critical impact of batch
  effects in high-throughput data. Nat Rev Genet, 11:733-739.
  doi:10.1038/nrg2825.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")

DATA_DIR = Path(__file__).parent.parent / "data"

INTERNEURON_SUBTYPES = {
    "PVALB": "PV+",
    "SST": "SST+",
    "NPY": "SST+",
    "PENK": "SST+",
    "VIP": "VIP+",
    "CALB2": "VIP+",
    "CALB1": "CB+",
    "GAD1": "Pan-GABA",
    "GAD2": "Pan-GABA",
    "CAMK2A": "Excitatory",
    "SLC17A7": "Excitatory",
    "CUX2": "Excitatory",
    "RORB": "Excitatory",
    "MBP": "Oligodendrocyte",
    "CNP": "Oligodendrocyte",
    "DLG4": "Synaptic",
    "VAMP1": "Synaptic",
    "AKT1": "Risk/Signaling",
    "COMT": "Risk/Signaling",
    "DISC1": "Risk/Signaling",
    "BDNF": "Risk/Signaling",
    "NR4A2": "Risk/Signaling",
    "RGS4": "Risk/Signaling",
    "RELN": "Risk/Signaling",
    "KCNH2": "Risk/Signaling",
    "CNR1": "Risk/Signaling",
    "PPP1R1B": "Risk/Signaling",
    "ERBB4": "Risk/Signaling",
    "ERBB3": "Risk/Signaling",
}

DEFAULT_SUBTYPE    = "SCZ candidate"
GABAERGIC_SUBTYPES = {"PV+", "SST+", "VIP+", "CB+", "Pan-GABA"}


def load_raw():
    path = DATA_DIR / "gse53987_dlpfc_raw.csv"
    if not path.exists():
        raise FileNotFoundError(str(path) + " not found. Run 01_fetch_geo.py first.")
    df = pd.read_csv(path)
    print("Loaded: " + str(df.shape[0]) + " records")
    return df


def filter_scz_control(df):
    """
    Keep only Schizophrenia and Control donors.
    GSE53987 contains four diagnostic groups: Schizophrenia, Control,
    Bipolar disorder, and Major depressive disorder. Including the non-SCZ
    psychiatric groups would confound the primary SCZ vs Control comparison
    and reduce statistical power for differential expression analysis.
    """
    before = df["donor_id"].nunique()
    df     = df[df["diagnosis"].isin(["Schizophrenia", "Control"])].copy()
    after  = df["donor_id"].nunique()
    n_scz  = df[df["diagnosis"] == "Schizophrenia"]["donor_id"].nunique()
    n_ctrl = df[df["diagnosis"] == "Control"]["donor_id"].nunique()

    print("\nFiltered to SCZ + Control: " + str(before) + " -> " + str(after) + " donors")
    print("  SCZ: " + str(n_scz) + "   Control: " + str(n_ctrl))
    return df


def zscore_within_gene(x):
    """
    Compute z-score for a single gene's expression values across donors.
    Used inside a groupby-transform call in transform().
    The 1e-9 term prevents division by zero for genes with zero variance.
    """
    return (x - x.mean()) / (x.std() + 1e-9)


def transform(df):
    """Log2 transform and flag outliers (> 3 SD within gene)."""

    # Search for the first valid expression column.
    # A flag stops the search once we find one with actual data.
    candidate_cols = ["cell_density", "staining_intensity", "expression_density"]
    primary        = None
    found_primary  = False

    for col in candidate_cols:
        if col in df.columns and df[col].notna().sum() > 0 and not found_primary:
            primary       = col
            found_primary = True

    if primary is None:
        raise ValueError("No valid expression column found.")

    df = df.dropna(subset=[primary]).copy()

    # Log2 transform: clip to 0 first so we never take log2 of a negative number.
    # Adding 1 before log2 (log2(x+1)) avoids log2(0) = -infinity for zero values.
    df["log2_expr"] = np.log2(df[primary].clip(lower=0) + 1)

    # Compute a z-score per gene across donors to identify extreme outliers.
    # groupby("gene").transform applies zscore_within_gene to each gene's
    # values separately, returning a Series with the same index as df.
    z = df.groupby("gene")["log2_expr"].transform(zscore_within_gene)

    df["is_outlier"] = z.abs() > 3

    print("\nTransformed: " + str(len(df)) + " records  |  outliers flagged: "
          + str(df["is_outlier"].sum()))
    return df


def build_pivot(df):
    """Aggregate to donor x gene matrix (mean log2 expression per donor-gene pair)."""

    # Average across any replicate probes for the same donor-gene pair.
    agg   = df.groupby(["donor_id", "gene"])["log2_expr"].mean().reset_index()
    pivot = agg.pivot(index="donor_id", columns="gene", values="log2_expr")

    # Fill any missing gene-donor combinations with the gene's mean across donors.
    pivot = pivot.fillna(pivot.mean())

    # Build a metadata DataFrame keeping only columns that actually exist.
    meta_cols = []
    for col in ["donor_id", "diagnosis", "age_days", "sex", "hemisphere"]:
        if col in df.columns:
            meta_cols.append(col)

    meta = (
        df[meta_cols].drop_duplicates("donor_id")
                     .set_index("donor_id")
                     .reindex(pivot.index)
    )

    return pivot, meta


def correct_batch(pivot, meta):
    """
    Linear batch correction by regressing out donor processing cohort.
    Batch is inferred from the GSM ID prefix:
      GSM1304xxx -> Batch A
      GSM1305xxx -> Batch B
    For each gene: fit expr ~ C(batch), take residuals + intercept.
    This removes the systematic between-batch expression shift while
    preserving within-batch biological variation.
    """
    try:
        from statsmodels.formula.api import ols
    except ImportError:
        print("  statsmodels not installed - skipping batch correction")
        print("  Run: pip install statsmodels")
        return pivot

    print("\nDetecting batches from donor ID prefix...")

    # Assign each donor to a batch based on their GSM ID prefix.
    # GSM1304xxx = Batch A, GSM1305xxx = Batch B.
    batch_labels = []
    for d in pivot.index:
        if d.startswith("GSM1305"):
            batch_labels.append("B")
        else:
            batch_labels.append("A")

    batch  = pd.Series(batch_labels, index=pivot.index, name="batch")
    counts = batch.value_counts()

    print("  Batch A (GSM1304): " + str(counts.get("A", 0)) + " donors")
    print("  Batch B (GSM1305): " + str(counts.get("B", 0)) + " donors")

    # If all donors are in the same batch, regression has nothing to correct.
    if counts.nunique() == 1 or len(counts) < 2:
        print("  Only one batch detected - skipping correction")
        return pivot

    corrected   = pivot.copy()
    n_corrected = 0

    for gene in pivot.columns:
        tmp = pd.DataFrame({"expr": pivot[gene], "batch": batch}).dropna()

        # Only attempt correction if both batches have data for this gene.
        if tmp["batch"].nunique() >= 2:
            try:
                # Fit a linear model: expression ~ batch label.
                # The residuals capture expression variation not explained by batch.
                # Adding back the intercept re-centers residuals at the grand mean.
                model     = ols("expr ~ C(batch)", data=tmp).fit()
                residuals = model.resid
                intercept = model.params["Intercept"]
                corrected.loc[tmp.index, gene] = residuals + intercept
                n_corrected = n_corrected + 1
            except Exception:
                pass  # if the model fails for a gene, leave it uncorrected

    print("  Batch corrected: " + str(n_corrected) + " genes")
    return corrected


def zscore(pivot):
    """Z-score normalize: each gene has mean=0, std=1 across donors."""
    scaler = StandardScaler()
    return pd.DataFrame(
        scaler.fit_transform(pivot),
        index=pivot.index,
        columns=pivot.columns
    )


def build_gene_annotations(genes):
    """Annotate each gene with interneuron subtype label."""

    subtypes      = []
    is_gabaergic  = []

    for gene in genes:
        subtype = INTERNEURON_SUBTYPES.get(gene, DEFAULT_SUBTYPE)
        subtypes.append(subtype)
        is_gabaergic.append(subtype in GABAERGIC_SUBTYPES)

    return pd.DataFrame({
        "gene":         genes,
        "subtype":      subtypes,
        "is_gabaergic": is_gabaergic,
    })


def save(expr_z, expr_log, meta, gene_annot):
    expr_z.to_csv(DATA_DIR / "expression_matrix.csv")
    expr_log.to_csv(DATA_DIR / "expression_matrix_raw.csv")
    meta.to_csv(DATA_DIR / "donor_metadata.csv")
    gene_annot.to_csv(DATA_DIR / "gene_annotations.csv", index=False)

    print("\nSaved to " + str(DATA_DIR) + "/")
    print("  expression_matrix.csv      (z-scored)")
    print("  expression_matrix_raw.csv  (log2)")
    print("  donor_metadata.csv")
    print("  gene_annotations.csv")


def main():
    print("=" * 60)
    print(" Step 2: Preprocessing")
    print("=" * 60)

    df = load_raw()
    df = filter_scz_control(df)
    df = transform(df)

    pivot, meta = build_pivot(df)
    print("\nExpression matrix shape: " + str(pivot.shape))

    # correct_batch() detects batch from GSM ID prefix and applies linear
    # regression correction. With DLPFC-only samples (n=34), all donors fall
    # in GSM1304xxx: one batch detected, correction is skipped automatically.
    # Function retained for use with multi-region data.
    pivot_bc = correct_batch(pivot, meta)

    expr_z    = zscore(pivot_bc)
    expr_log  = pivot_bc  # log2-transformed; batch correction no-ops on single-cohort DLPFC data

    gene_annot = build_gene_annotations(list(expr_z.columns))

    save(expr_z, expr_log, meta, gene_annot)

    if "diagnosis" in meta.columns:
        counts = meta["diagnosis"].value_counts()
        for label in counts.index:
            print("  " + str(label) + ": " + str(counts[label]) + " donors")

    print("\nDone. Run 03_cluster_analysis.py next.")


if __name__ == "__main__":
    main()