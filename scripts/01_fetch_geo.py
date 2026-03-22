"""
01_fetch_geo.py
===============
Downloads postmortem DLPFC gene expression data from NCBI GEO.

Dataset: GSE53987
  Tissue:    Postmortem prefrontal cortex (Brodmann Area 46), striatum, hippocampus
  Platform:  GPL570 — Affymetrix Human Genome U133 Plus 2.0 Array
  Groups:    Schizophrenia, bipolar disorder, major depressive disorder, controls
  Reference: Lanz TA, et al. Translational Psychiatry (2019). PMID: 31123247
  URL:       https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53987

This script downloads via HTTPS (not FTP) to avoid connection drop issues,
and parses the SOFT file directly rather than relying on GEOparse's downloader.
GEOparse was attempted first but its FTP downloader produced incomplete files
due to connection timeouts; HTTPS via requests is more reliable.

This analysis uses only the prefrontal cortex schizophrenia and control samples.
Bipolar disorder and major depressive disorder samples are excluded in
02_preprocess.py.

Gene panel note
---------------
TARGET_GENES below lists the 58-gene panel from Guillozet-Bongaarts et al. (2014),
which this project was designed to replicate computationally. Not all genes are
present on the GPL570 platform — genes absent from the platform (including CHRNA7
and PRODH) will be reported as "Not in platform" at runtime and excluded from
analysis. 48 of the 58 target genes were found on the platform.

Usage
-----
  # Delete any incomplete cached files before first run:
  rm -rf data/geo_cache/

  python scripts/01_fetch_geo.py

Outputs
-------
  data/allen_scz_raw.csv        Long-format expression (primary input for pipeline)
  data/geo_expression_wide.csv  Wide matrix (genes x donors)
  data/geo_metadata.csv         Sample metadata

References
----------
Lanz TA, Reinhart V, Sheehan MJ, Rizzo SJS, et al. (2019). Postmortem
  transcriptional profiling reveals widespread increase in inflammation in
  schizophrenia. Translational Psychiatry, 9(1):151.
  doi:10.1038/s41398-019-0494-6. PMID: 31123247. GEO: GSE53987.

Guillozet-Bongaarts AL, Hyde TM, Dalley RA, et al. (2014). Altered gene
  expression in the dorsolateral prefrontal cortex of individuals with
  schizophrenia. Molecular Psychiatry, 19(4):478-485.
  doi:10.1038/mp.2013.30. PMID: 23528911.
"""

import sys
import gzip
import shutil
import requests
import pandas as pd
import numpy as np
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
CACHE_DIR = DATA_DIR / "geo_cache"
DATA_DIR.mkdir(exist_ok=True)
CACHE_DIR.mkdir(exist_ok=True)

GEO_ACCESSION = "GSE53987"
SOFT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/"
    "GSE53987/soft/GSE53987_family.soft.gz"
)
SOFT_PATH = CACHE_DIR / "GSE53987_family.soft.gz"
SOFT_TEXT  = CACHE_DIR / "GSE53987_family.soft"

# Target genes from the 58-gene panel of Guillozet-Bongaarts et al. (2014).
# Note: some genes (e.g. CHRNA7, PRODH) are absent from the GPL570 platform
# and will be reported as "Not in platform" at runtime. 48 of 58 were found.
TARGET_GENES = [
    "GAD1","GAD2","PVALB","SST","NPY","VIP","CALB1","CALB2","RELN",
    "PENK","DISC1","BDNF","COMT","NR4A2","RGS4","MBP","DLG4","AKT1",
    "ERBB4","ERBB3","RORB","CUX2","SLC17A7","CHRNA7","CNR1","PRODH",
    "KCNH2","CNP","CAMK2A","PPP1R1B","NEFH","PCP4","GAP43","SLC1A2",
    "SLC6A1","TAC1","NDEL1","NOS1AP","SNCG","SYNPR","VAMP1","CLDN5",
    "NTNG2","NDNF","CARTPT","FEZ1","RASGRF2","CTGF","MFGE8","BLOC1S6",
]


# ── Download via HTTPS ────────────────────────────────────────────────────────
def download_soft():
    """Download the SOFT file via HTTPS with a progress bar."""
    if SOFT_TEXT.exists():
        print(f"Using cached SOFT file: {SOFT_TEXT}")
        return

    if SOFT_PATH.exists():
        print(f"Deleting incomplete cached gz: {SOFT_PATH}")
        SOFT_PATH.unlink()

    print(f"Downloading {GEO_ACCESSION} from NCBI GEO (HTTPS)...")
    print(f"URL: {SOFT_URL}\n")

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
            "AppleWebKit/537.36 Chrome/120.0.0.0 Safari/537.36"
        )
    }

    with requests.get(SOFT_URL, headers=headers, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        downloaded = 0
        with open(SOFT_PATH, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 256):
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded / total * 100
                    mb  = downloaded / 1e6
                    print(f"\r  {mb:.1f} MB / {total/1e6:.1f} MB  ({pct:.0f}%)",
                          end="", flush=True)
    print(f"\n  Download complete: {SOFT_PATH.stat().st_size / 1e6:.1f} MB")

    print("  Decompressing...")
    with gzip.open(SOFT_PATH, "rb") as gz, open(SOFT_TEXT, "wb") as out:
        shutil.copyfileobj(gz, out)
    print(f"  Decompressed: {SOFT_TEXT.stat().st_size / 1e6:.1f} MB")


# ── Parse SOFT file ───────────────────────────────────────────────────────────
def parse_soft() -> tuple:
    """
    Parse the GEO SOFT file into:
      - metadata dict   {sample_id -> {key -> value}}
      - expression dict {sample_id -> {probe_id -> value}}
      - platform dict   {probe_id  -> gene_symbol}
    """
    print("\nParsing SOFT file...")

    metadata   = {}   # sample_id -> {key: value}
    expression = {}   # sample_id -> {probe: value}
    platform   = {}   # probe_id  -> gene_symbol

    current_sample    = None
    current_platform  = None
    in_sample_table   = False
    in_platform_table = False
    platform_header   = []
    sym_col_idx       = None

    sym_candidates = [
        "Gene Symbol","GENE_SYMBOL","gene_symbol","Symbol",
        "SYMBOL","Gene_Symbol","GENE SYMBOL",
    ]

    with open(SOFT_TEXT, "r", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")

            # ── Section markers ────────────────────────────────────────────
            if line.startswith("^SAMPLE"):
                current_sample = line.split("=", 1)[1].strip()
                metadata[current_sample] = {}
                in_sample_table = False
                continue

            if line.startswith("^PLATFORM"):
                current_platform = line.split("=", 1)[1].strip()
                in_platform_table = False
                continue

            # ── Platform table (probe → gene mapping) ─────────────────────
            if line == "!platform_table_begin":
                in_platform_table = True
                platform_header = []
                sym_col_idx = None
                continue

            if line == "!platform_table_end":
                in_platform_table = False
                continue

            if in_platform_table:
                parts = line.split("\t")
                if not platform_header:
                    platform_header = parts
                    sym_col_idx = next(
                        (i for i, c in enumerate(parts) if c in sym_candidates),
                        None
                    )
                    if sym_col_idx is None:
                        print(f"  Platform cols: {parts[:8]}")
                    continue
                if sym_col_idx is not None and len(parts) > sym_col_idx:
                    probe_id = parts[0].strip()
                    gene_sym = parts[sym_col_idx].strip()
                    if gene_sym and probe_id:
                        # Handle "GENE1 /// GENE2" — take first
                        gene_sym = gene_sym.split("///")[0].strip().upper()
                        if gene_sym:
                            platform[probe_id] = gene_sym
                continue

            # ── Sample table (expression values) ──────────────────────────
            if line == "!sample_table_begin":
                in_sample_table = True
                if current_sample not in expression:
                    expression[current_sample] = {}
                continue

            if line == "!sample_table_end":
                in_sample_table = False
                continue

            if in_sample_table and current_sample:
                parts = line.split("\t")
                if len(parts) >= 2 and parts[0] != "ID_REF":
                    try:
                        expression[current_sample][parts[0]] = float(parts[1])
                    except (ValueError, IndexError):
                        pass
                continue

            # ── Sample metadata ────────────────────────────────────────────
            if (current_sample and
                    line.startswith("!Sample_") and
                    not in_sample_table):
                if "=" in line:
                    k, v = line.split("=", 1)
                    k = k.replace("!Sample_", "").strip().lower()
                    v = v.strip()
                    # Append repeated keys (e.g. characteristics_ch1)
                    if k in metadata[current_sample]:
                        existing = metadata[current_sample][k]
                        if isinstance(existing, list):
                            existing.append(v)
                        else:
                            metadata[current_sample][k] = [existing, v]
                    else:
                        metadata[current_sample][k] = v

    n_samples  = len(expression)
    n_probes   = len(platform)
    probe_vals = sum(len(v) for v in expression.values())
    print(f"  Samples: {n_samples}")
    print(f"  Platform probes with gene symbols: {n_probes}")
    print(f"  Total probe-value pairs: {probe_vals:,}")

    return metadata, expression, platform


# ── Build expression matrix ───────────────────────────────────────────────────
def build_matrix(metadata: dict, expression: dict, platform: dict) -> tuple:
    """
    Build gene-level expression matrix and sample metadata DataFrame.

    For genes with multiple probes, the probe with highest mean expression
    across all samples is selected (a standard approach for Affymetrix data).
    """
    # ── Metadata DataFrame ────────────────────────────────────────────────
    meta_rows = []
    for sample_id, m in metadata.items():
        row = {"sample_id": sample_id}

        # Flatten characteristics_ch1 list
        chars = m.get("characteristics_ch1", [])
        if isinstance(chars, str):
            chars = [chars]
        for ch in chars:
            if ":" in ch:
                k, v = ch.split(":", 1)
                row[k.strip().lower().replace(" ", "_")] = v.strip()

        # Basic fields
        for key in ["title","geo_accession","source_name_ch1",
                    "organism_ch1","molecule_ch1"]:
            val = m.get(key, "")
            if isinstance(val, list):
                val = val[0]
            row[key] = val

        meta_rows.append(row)

    meta = pd.DataFrame(meta_rows).set_index("sample_id")

    # Find diagnosis column
    meta["diagnosis"] = "Unknown"
    for col in meta.columns:
        if col in ["sample_id","diagnosis"]:
            continue
        try:
            vals = meta[col].dropna().astype(str).str.lower()
            if vals.str.contains(r"schizophreni|control", regex=True, na=False).any():
                meta["diagnosis"] = (
                    meta[col].astype(str).str.strip()
                    .str.replace(r"(?i).*schizophreni.*", "Schizophrenia", regex=True)
                    .str.replace(r"(?i).*control.*",      "Control",       regex=True)
                )
                print(f"\n  Diagnosis from column '{col}':")
                print(f"  {meta['diagnosis'].value_counts().to_dict()}")
                break
        except Exception:
            continue

    # Age
    for col in meta.columns:
        if "age" in col.lower() and col != "diagnosis":
            try:
                meta["age_years"] = pd.to_numeric(meta[col], errors="coerce")
                meta["age_days"]  = (meta["age_years"] * 365.25).round()
                print(f"  Age from column '{col}'")
                break
            except Exception:
                pass

    # Sex
    for col in ["sex","gender","Sex","Gender"]:
        if col in meta.columns:
            meta["sex"] = meta[col].astype(str).str.strip()
            break

    print(f"  Total samples in metadata: {len(meta)}")

    # ── Expression matrix ─────────────────────────────────────────────────
    # Map probe IDs to gene symbols, keep only target genes
    target_set = set(TARGET_GENES)
    target_probes = {
        probe: gene
        for probe, gene in platform.items()
        if gene in target_set
    }

    print(f"\nBuilding gene expression matrix...")
    print(f"  Target probes mapped: {len(target_probes)}")

    # Collect expression values per gene across all samples
    gene_probe_vals = {}  # gene -> {probe -> {sample_id -> value}}
    for sample_id, probe_vals in expression.items():
        if sample_id not in meta.index:
            continue
        for probe, val in probe_vals.items():
            gene = target_probes.get(probe)
            if gene is None:
                continue
            if gene not in gene_probe_vals:
                gene_probe_vals[gene] = {}
            if probe not in gene_probe_vals[gene]:
                gene_probe_vals[gene][probe] = {}
            gene_probe_vals[gene][probe][sample_id] = val

    # For each gene, pick the probe with highest mean expression
    sample_ids = list(meta.index)
    gene_expr  = {}  # gene -> pd.Series indexed by sample_id

    for gene, probes in gene_probe_vals.items():
        best_probe = None
        best_mean  = -np.inf
        for probe, sample_dict in probes.items():
            vals  = [sample_dict.get(s, np.nan) for s in sample_ids]
            nmean = np.nanmean(vals) if any(not np.isnan(v) for v in vals) else -np.inf
            if nmean > best_mean:
                best_mean  = nmean
                best_probe = probe
        vals = [gene_probe_vals[gene][best_probe].get(s, np.nan) for s in sample_ids]
        gene_expr[gene] = pd.Series(vals, index=sample_ids)

    expr_wide = pd.DataFrame(gene_expr).T   # genes x samples
    expr_wide.index.name = "gene"

    found   = sorted(set(TARGET_GENES) & set(expr_wide.index))
    missing = sorted(set(TARGET_GENES) - set(expr_wide.index))
    print(f"  Genes found: {len(found)} / {len(TARGET_GENES)}")
    if missing:
        print(f"  Not in platform: {missing}")
    if found:
        expr_wide = expr_wide.loc[found]

    return expr_wide, meta


# ── Reshape to long format ────────────────────────────────────────────────────
def to_long(expr: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Reshape wide expression matrix to long format for downstream pipeline."""
    print("\nReshaping to long format...")
    long = (
        expr.reset_index()
            .melt(id_vars="gene", var_name="donor_id", value_name="cell_density")
    )
    meta_reset = meta.reset_index().rename(columns={"sample_id": "donor_id"})
    keep = ["donor_id"] + [
        c for c in ["diagnosis","age_years","age_days","sex"]
        if c in meta_reset.columns
    ]
    long = long.merge(meta_reset[keep], on="donor_id", how="left")
    if "age_days" not in long.columns and "age_years" in long.columns:
        long["age_days"] = (long["age_years"] * 365.25).round()
    long["structure_name"]    = "dorsolateral prefrontal cortex"
    long["structure_acronym"] = "DLPFC"
    long["hemisphere"]        = "unknown"
    print(f"  Shape: {long.shape}")
    return long


# ── Save ──────────────────────────────────────────────────────────────────────
def save(long: pd.DataFrame, wide: pd.DataFrame, meta: pd.DataFrame):
    long.to_csv(DATA_DIR / "allen_scz_raw.csv",        index=False)
    wide.to_csv(DATA_DIR / "geo_expression_wide.csv")
    meta.to_csv(DATA_DIR / "geo_metadata.csv")

    print(f"\n{'='*60}")
    print("SUCCESS")
    print(f"{'='*60}")
    print(f"  data/allen_scz_raw.csv   ({len(long):,} records)")
    print(f"  Genes   : {long['gene'].nunique()}")
    print(f"  Donors  : {long['donor_id'].nunique()}")
    if "diagnosis" in long.columns:
        deduped = long.drop_duplicates("donor_id")
        for k, v in deduped["diagnosis"].value_counts().items():
            print(f"  {k}: {v} donors")
    print(f"\nDataset: {GEO_ACCESSION} — DLPFC SCZ vs Control (Lanz et al. 2019)")
    print(f"\nNext: python scripts/02_preprocess.py")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print(f" GEO Fetcher — {GEO_ACCESSION}")
    print(f" DLPFC Schizophrenia vs Control")
    print("=" * 60)

    download_soft()
    metadata, expression, platform = parse_soft()

    if not expression:
        print("ERROR: No expression data parsed from SOFT file.")
        sys.exit(1)

    wide, meta = build_matrix(metadata, expression, platform)

    if wide.empty:
        print("ERROR: Expression matrix is empty.")
        sys.exit(1)

    long = to_long(wide, meta)
    save(long, wide, meta)


if __name__ == "__main__":
    main()
