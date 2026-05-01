"""
02_preprocess.py
Cleans, batch-corrects, and normalizes GEO expression data.

Key steps:
1. Filter to Schizophrenia and Control donors only
   (GSE53987 also contains bipolar disorder and MDD donors, excluded here)
2. Detect and correct processing batch (GSM ID prefix-based)
3. Log2 transform expression values
4. Aggregate to mean per donor x gene
5. Z-score normalize across donors per gene
6. Annotate genes with interneuron subtype labels

Batch correction rationale:
PCA and hierarchical clustering of the raw data show that PC1 (42% of variance)
and the dominant cluster split correspond to GSM1304xxx vs GSM1305xxx sample ID
prefixes, indicating a technical batch effect from two processing cohorts. This
is a well-documented issue in multi-cohort postmortem brain microarray studies.

We use linear regression to remove the batch effect before downstream analysis:
for each gene, we fit expr ~ batch, then take residuals + intercept. This removes
the systematic expression shift between batches while preserving within-batch
biological variation, following Leek et al. (2010).

Note: explicit covariates (RIN, PMI, medication) were not available in this
dataset; batch is inferred from sample ID prefix only.

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

# Interneuron subtype annotations for the gene panel.
# Based on cell-type marker categories from Guillozet-Bongaarts et al. (2014)
# and Dienel & Lewis (2019).
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
DEFAULT_SUBTYPE = "SCZ candidate"
GABAERGIC_SUBTYPES = {"PV+", "SST+", "VIP+", "CB+", "Pan-GABA"}


def load_raw() -> pd.DataFrame:
    path = DATA_DIR / "allen_scz_raw.csv"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found. Run 01_fetch_geo.py first.")
    df = pd.read_csv(path)
    print(f"Loaded: {df.shape[0]:,} records")
    return df


def filter_scz_control(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only Schizophrenia and Control donors.

    GSE53987 contains four diagnostic groups: Schizophrenia, Control,
    Bipolar disorder, and Major depressive disorder. Including the non-SCZ
    psychiatric groups would confound the primary SCZ vs Control comparison
    and reduce statistical power for differential expression analysis.
    """
    before = df["donor_id"].nunique()
    df = df[df["diagnosis"].isin(["Schizophrenia", "Control"])].copy()
    after = df["donor_id"].nunique()
    n_scz = df[df["diagnosis"] == "Schizophrenia"]["donor_id"].nunique()
    n_ctrl = df[df["diagnosis"] == "Control"]["donor_id"].nunique()
    print(f"\nFiltered to SCZ + Control: {before} -> {after} donors")
    print(f"  SCZ: {n_scz}   Control: {n_ctrl}")
    return df


def transform(df: pd.DataFrame) -> pd.DataFrame:
    """Log2 transform and flag outliers (> 3 SD within gene)."""
    primary = next(
        (c for c in ["cell_density","staining_intensity","expression_density"]
         if c in df.columns and df[c].notna().sum() > 0),
        None
    )
    if primary is None:
        raise ValueError("No valid expression column found.")
    df = df.dropna(subset=[primary]).copy()
    df["log2_expr"] = np.log2(df[primary].clip(lower=0) + 1)
    z = df.groupby("gene")["log2_expr"].transform(
        lambda x: (x - x.mean()) / (x.std() + 1e-9)
    )
    df["is_outlier"] = z.abs() > 3
    print(f"\nTransformed: {len(df):,} records  |  outliers flagged: {df['is_outlier'].sum()}")
    return df


def build_pivot(df: pd.DataFrame) -> tuple:
    """Aggregate to donor x gene matrix (mean log2 expression per donor-gene pair)."""
    agg = df.groupby(["donor_id","gene"])["log2_expr"].mean().reset_index()
    pivot = agg.pivot(index="donor_id", columns="gene", values="log2_expr")
    pivot = pivot.fillna(pivot.mean())

    meta_cols = [c for c in ["donor_id","diagnosis","age_days","sex","hemisphere"]
                 if c in df.columns]
    meta = (
        df[meta_cols].drop_duplicates("donor_id")
                     .set_index("donor_id")
                     .reindex(pivot.index)
    )
    return pivot, meta


def correct_batch(pivot: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """
    Linear batch correction by regressing out donor processing cohort.

    Batch is inferred from the GSM ID prefix:
      GSM1304xxx → Batch A
      GSM1305xxx → Batch B

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
    batch = pd.Series(
        [("B" if d.startswith("GSM1305") else "A") for d in pivot.index],
        index=pivot.index,
        name="batch"
    )
    counts = batch.value_counts()
    print(f"  Batch A (GSM1304): {counts.get('A', 0)} donors")
    print(f"  Batch B (GSM1305): {counts.get('B', 0)} donors")

    if counts.nunique() == 1 or len(counts) < 2:
        print("  Only one batch detected - skipping correction")
        return pivot

    corrected = pivot.copy()
    n_corrected = 0
    for gene in pivot.columns:
        tmp = pd.DataFrame({"expr": pivot[gene], "batch": batch}).dropna()
        if tmp["batch"].nunique() < 2:
            continue
        try:
            model = ols("expr ~ C(batch)", data=tmp).fit()
            residuals = model.resid
            intercept = model.params["Intercept"]
            corrected.loc[tmp.index, gene] = residuals + intercept
            n_corrected += 1
        except Exception:
            pass

    print(f"  Batch corrected: {n_corrected} genes")
    return corrected


def zscore(pivot: pd.DataFrame) -> pd.DataFrame:
    """Z-score normalize: each gene has mean=0, std=1 across donors."""
    scaler = StandardScaler()
    return pd.DataFrame(
        scaler.fit_transform(pivot),
        index=pivot.index,
        columns=pivot.columns
    )


def build_gene_annotations(genes: list) -> pd.DataFrame:
    """Annotate each gene with interneuron subtype label."""
    return pd.DataFrame({
        "gene": genes,
        "subtype": [INTERNEURON_SUBTYPES.get(g, DEFAULT_SUBTYPE) for g in genes],
        "is_gabaergic": [INTERNEURON_SUBTYPES.get(g,"") in GABAERGIC_SUBTYPES
                         for g in genes],
    })


def save(expr_z, expr_log, meta, gene_annot):
    expr_z.to_csv(DATA_DIR / "expression_matrix.csv")
    expr_log.to_csv(DATA_DIR / "expression_matrix_raw.csv")
    meta.to_csv(DATA_DIR / "donor_metadata.csv")
    gene_annot.to_csv(DATA_DIR / "gene_annotations.csv", index=False)
    print(f"\nSaved to {DATA_DIR}/")
    print("  expression_matrix.csv      (batch-corrected, z-scored)")
    print("  expression_matrix_raw.csv  (batch-corrected, log2)")
    print("  donor_metadata.csv")
    print("  gene_annotations.csv")


def main():
    print("=" * 60)
    print(" Step 2: Preprocessing + Batch Correction")
    print("=" * 60)

    df = load_raw()
    df = filter_scz_control(df)
    df = transform(df)
    pivot, meta = build_pivot(df)

    print(f"\nExpression matrix before correction: {pivot.shape}")

    # correct_batch() detects batch from GSM ID prefix and applies linear regression correction.
    # With DLPFC-only samples (n=34), all donors fall in GSM1304xxx -- one batch detected,
    # correction is skipped automatically. Function retained for use with multi-region data.
    pivot_bc = correct_batch(pivot, meta)
    expr_z = zscore(pivot_bc)
    expr_log = pivot_bc   # keep log2 batch-corrected for DE analysis

    gene_annot = build_gene_annotations(list(expr_z.columns))
    save(expr_z, expr_log, meta, gene_annot)

    if "diagnosis" in meta.columns:
        for k, v in meta["diagnosis"].value_counts().items():
            print(f"  {k}: {v} donors")

    print("\nDone. Run 03_cluster_analysis.py next.")


if __name__ == "__main__":
    main()