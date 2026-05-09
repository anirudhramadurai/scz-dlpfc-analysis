"""
01_fetch_geo.py
Downloads postmortem DLPFC gene expression data from NCBI GEO.

Dataset: GSE53987
  Tissue:    Postmortem prefrontal cortex (Brodmann Area 46), striatum, hippocampus
  Platform:  GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array
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

Gene panel note:
TARGET_GENES below lists 50 genes drawn from the 58-gene panel in Guillozet-Bongaarts et al. (2014).
Of those 50, CHRNA7 and PRODH are absent from the GPL570 platform, will be reported as
"Not in platform" at runtime, and excluded from analysis, leaving 48 genes used.
The remaining 8 genes from the source panel were not carried forward into TARGET_GENES.

Usage:
  # Delete any incomplete cached files before first run:
  rm -rf data/geo_cache/

  python scripts/01_fetch_geo.py

Outputs:
  data/allen_scz_raw.csv        Long-format expression (primary input for pipeline)
  data/geo_expression_wide.csv  Wide matrix (genes x donors)
  data/geo_metadata.csv         Sample metadata

References:
Lanz TA, Reinhart V, Sheehan MJ, Rizzo SJS, et al. (2019). Postmortem
  transcriptional profiling reveals widespread increase in inflammation in
  schizophrenia. Translational Psychiatry, 9(1):151.
  doi:10.1038/s41398-019-0492-8. PMID: 31123247. GEO: GSE53987.

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

DATA_DIR  = Path(__file__).parent.parent / "data"
CACHE_DIR = DATA_DIR / "geo_cache"
DATA_DIR.mkdir(exist_ok=True)
CACHE_DIR.mkdir(exist_ok=True)

GEO_ACCESSION = "GSE53987"

SOFT_URL  = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/"
    "GSE53987/soft/GSE53987_family.soft.gz"
)
SOFT_PATH = CACHE_DIR / "GSE53987_family.soft.gz"
SOFT_TEXT = CACHE_DIR / "GSE53987_family.soft"

TARGET_GENES = [
    "GAD1","GAD2","PVALB","SST","NPY","VIP","CALB1","CALB2","RELN",
    "PENK","DISC1","BDNF","COMT","NR4A2","RGS4","MBP","DLG4","AKT1",
    "ERBB4","ERBB3","RORB","CUX2","SLC17A7","CHRNA7","CNR1","PRODH",
    "KCNH2","CNP","CAMK2A","PPP1R1B","NEFH","PCP4","GAP43","SLC1A2",
    "SLC6A1","TAC1","NDEL1","NOS1AP","SNCG","SYNPR","VAMP1","CLDN5",
    "NTNG2","NDNF","CARTPT","FEZ1","RASGRF2","CTGF","MFGE8","BLOC1S6",
]


def download_soft():
    """Download the SOFT file via HTTPS with a progress bar."""

    if SOFT_TEXT.exists():
        print("Using cached SOFT file: " + str(SOFT_TEXT))
    else:
        if SOFT_PATH.exists():
            print("Deleting incomplete cached gz: " + str(SOFT_PATH))
            SOFT_PATH.unlink()

        print("Downloading " + GEO_ACCESSION + " from NCBI GEO (HTTPS)...")
        print("URL: " + SOFT_URL + "\n")

        headers = {
            "User-Agent": (
                "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                "AppleWebKit/537.36 Chrome/120.0.0.0 Safari/537.36"
            )
        }

        with requests.get(SOFT_URL, headers=headers, stream=True, timeout=120) as r:
            r.raise_for_status()
            total      = int(r.headers.get("content-length", 0))
            downloaded = 0

            with open(SOFT_PATH, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 256):
                    f.write(chunk)
                    downloaded = downloaded + len(chunk)
                    if total > 0:
                        pct = downloaded / total * 100
                        mb  = downloaded / 1e6
                        print("\r  " + str(round(mb, 1)) + " MB / "
                              + str(round(total / 1e6, 1)) + " MB  ("
                              + str(round(pct)) + "%)",
                              end="", flush=True)

        size_mb = round(SOFT_PATH.stat().st_size / 1e6, 1)
        print("\n  Download complete: " + str(size_mb) + " MB")

        print("  Decompressing...")
        with gzip.open(SOFT_PATH, "rb") as gz, open(SOFT_TEXT, "wb") as out:
            shutil.copyfileobj(gz, out)

        decompressed_mb = round(SOFT_TEXT.stat().st_size / 1e6, 1)
        print("  Decompressed: " + str(decompressed_mb) + " MB")


def parse_soft():
    """
    Parse the GEO SOFT file into three dictionaries:
      metadata   -- sample_id -> {key -> value}
      expression -- sample_id -> {probe_id -> float}
      platform   -- probe_id  -> gene_symbol

    The SOFT file has sections separated by lines starting with ^.
    Inside each section, metadata lines start with ! and data tables
    are bracketed by *_table_begin and *_table_end markers.

    Each line belongs to exactly one case, so the main loop is written
    as a single if/elif chain -- only one branch runs per line.
    """
    print("\nParsing SOFT file...")

    metadata   = {}
    expression = {}
    platform   = {}

    current_sample    = None
    in_sample_table   = False
    in_platform_table = False
    platform_header   = []
    sym_col_idx       = None

    sym_candidates = [
        "Gene Symbol", "GENE_SYMBOL", "gene_symbol", "Symbol",
        "SYMBOL", "Gene_Symbol", "GENE SYMBOL",
    ]

    with open(SOFT_TEXT, "r", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")

            # Each line falls into exactly one category, handled by elif so
            # only one branch runs. This replaces the original if + continue pattern.

            if line.startswith("^SAMPLE"):
                # Start of a new sample section.
                current_sample = line.split("=", 1)[1].strip()
                metadata[current_sample] = {}
                in_sample_table = False

            elif line.startswith("^PLATFORM"):
                # Start of the probe annotation section.
                in_platform_table = False

            elif line == "!platform_table_begin":
                in_platform_table = True
                platform_header   = []
                sym_col_idx       = None

            elif line == "!platform_table_end":
                in_platform_table = False

            elif in_platform_table:
                parts = line.split("\t")

                if len(platform_header) == 0:
                    # This is the header row. Find which column holds gene symbols.
                    platform_header = parts

                    # Use a flag so we stop updating sym_col_idx after the first match.
                    found_sym_col = False
                    for i in range(len(parts)):
                        if parts[i] in sym_candidates and not found_sym_col:
                            sym_col_idx   = i
                            found_sym_col = True

                    if sym_col_idx is None:
                        print("  Platform cols: " + str(parts[:8]))

                else:
                    # This is a data row: map probe_id -> gene_symbol.
                    if sym_col_idx is not None and len(parts) > sym_col_idx:
                        probe_id = parts[0].strip()
                        gene_sym = parts[sym_col_idx].strip()

                        if gene_sym != "" and probe_id != "":
                            # Some probes list multiple genes as "GENE1 /// GENE2".
                            # Take only the first gene for a 1-to-1 mapping.
                            gene_sym = gene_sym.split("///")[0].strip().upper()
                            if gene_sym != "":
                                platform[probe_id] = gene_sym

            elif line == "!sample_table_begin":
                in_sample_table = True
                if current_sample not in expression:
                    expression[current_sample] = {}

            elif line == "!sample_table_end":
                in_sample_table = False

            elif in_sample_table and current_sample is not None:
                # Data row inside a sample's expression table.
                # parts[0] = probe ID, parts[1] = raw expression value.
                parts = line.split("\t")
                if len(parts) >= 2 and parts[0] != "ID_REF":
                    try:
                        expression[current_sample][parts[0]] = float(parts[1])
                    except (ValueError, IndexError):
                        pass

            elif (current_sample is not None
                  and line.startswith("!Sample_")
                  and not in_sample_table):
                # Metadata line, e.g.: !Sample_characteristics_ch1 = tissue: Pre-frontal cortex
                if "=" in line:
                    k, v = line.split("=", 1)
                    k = k.replace("!Sample_", "").strip().lower()
                    v = v.strip()

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
    probe_vals = 0
    for sample_id in expression:
        probe_vals = probe_vals + len(expression[sample_id])

    print("  Samples: " + str(n_samples))
    print("  Platform probes with gene symbols: " + str(n_probes))
    print("  Total probe-value pairs: " + str(probe_vals))

    return metadata, expression, platform


def build_matrix(metadata, expression, platform):
    """
    Build a gene-level expression matrix and sample metadata DataFrame.

    Steps:
    1. Flatten raw metadata dicts into a tidy DataFrame.
    2. Parse diagnosis, age, and sex from GEO characteristics fields.
    3. Filter to DLPFC samples only.
    4. For each target gene, find all its probes and keep the one with
       the highest mean expression (standard Affymetrix selection strategy).
    5. Return a genes-x-samples DataFrame and the filtered metadata DataFrame.
    """

    # --- Build metadata DataFrame ---

    meta_rows = []
    for sample_id in metadata:
        m   = metadata[sample_id]
        row = {}
        row["sample_id"] = sample_id

        chars = m.get("characteristics_ch1", [])
        if isinstance(chars, str):
            chars = [chars]
        for ch in chars:
            if ":" in ch:
                k, v     = ch.split(":", 1)
                col_name = k.strip().lower().replace(" ", "_")
                row[col_name] = v.strip()

        for key in ["title", "geo_accession", "source_name_ch1",
                    "organism_ch1", "molecule_ch1"]:
            val = m.get(key, "")
            if isinstance(val, list):
                val = val[0]
            row[key] = val

        meta_rows.append(row)

    meta = pd.DataFrame(meta_rows).set_index("sample_id")

    # Extract diagnosis
    # Scan columns for one that contains both "schizophrenia" and "control".
    # found_diagnosis is a flag so we stop searching once we find the right column.
    meta["diagnosis"] = "Unknown"
    found_diagnosis   = False

    for col in meta.columns:
        if col != "sample_id" and col != "diagnosis" and not found_diagnosis:
            try:
                vals        = meta[col].dropna().astype(str).str.lower()
                has_scz     = vals.str.contains("schizophreni", regex=False, na=False).any()
                has_control = vals.str.contains("control",      regex=False, na=False).any()

                if has_scz and has_control:
                    meta["diagnosis"] = (
                        meta[col].astype(str).str.strip()
                        .str.replace(r"(?i).*schizophreni.*", "Schizophrenia", regex=True)
                        .str.replace(r"(?i).*control.*",      "Control",       regex=True)
                    )
                    print("\n  Diagnosis from column '" + col + "':")
                    print("  " + str(meta["diagnosis"].value_counts().to_dict()))
                    found_diagnosis = True

            except Exception:
                pass

    # Extract age
    found_age = False

    for col in meta.columns:
        if "age" in col.lower() and col != "diagnosis" and not found_age:
            try:
                meta["age_years"] = pd.to_numeric(meta[col], errors="coerce")
                meta["age_days"]  = (meta["age_years"] * 365.25).round()
                print("  Age from column '" + col + "'")
                found_age = True
            except Exception:
                pass

    # Extract sex
    found_sex    = False
    sex_cols     = ["sex", "gender", "Sex", "Gender"]

    for col in sex_cols:
        if col in meta.columns and not found_sex:
            meta["sex"] = meta[col].astype(str).str.strip()
            found_sex   = True

    print("  Total samples in metadata: " + str(len(meta)))

    # Filter to DLPFC only
    if "tissue" in meta.columns:
        dlpfc_mask = meta["tissue"].astype(str).str.contains(
            "Pre-frontal cortex|prefrontal|BA46", case=False, regex=True, na=False
        )
        meta = meta[dlpfc_mask].copy()
        print("  After DLPFC filter: " + str(len(meta)) + " samples")
        print("  Tissue values retained: " + str(meta["tissue"].unique().tolist()))
    else:
        print("  WARNING: 'tissue' column not found: no region filter applied")

    # Build a restricted probe->gene lookup for target genes only
    target_set    = set(TARGET_GENES)
    target_probes = {}
    for probe in platform:
        gene = platform[probe]
        if gene in target_set:
            target_probes[probe] = gene

    print("\nBuilding gene expression matrix...")
    print("  Target probes mapped: " + str(len(target_probes)))

    # Collect expression values grouped by gene and probe
    # Structure: gene -> probe -> sample_id -> float
    gene_probe_vals = {}

    for sample_id in expression:
        if sample_id in meta.index:
            probe_vals = expression[sample_id]
            for probe in probe_vals:
                gene = target_probes.get(probe)
                if gene is not None:
                    if gene not in gene_probe_vals:
                        gene_probe_vals[gene] = {}
                    if probe not in gene_probe_vals[gene]:
                        gene_probe_vals[gene][probe] = {}
                    gene_probe_vals[gene][probe][sample_id] = probe_vals[probe]

    # Select the best probe per gene
    # For each gene, pick the probe with the highest mean expression across
    # all DLPFC samples. Higher mean = better hybridization = more reliable signal.
    sample_ids = list(meta.index)
    gene_expr  = {}

    for gene in gene_probe_vals:
        probes     = gene_probe_vals[gene]
        best_probe = None
        best_mean  = -np.inf

        for probe in probes:
            sample_dict = probes[probe]

            vals = []
            for s in sample_ids:
                if s in sample_dict:
                    vals.append(sample_dict[s])
                else:
                    vals.append(np.nan)

            has_data = any(not np.isnan(v) for v in vals)
            if has_data:
                mean_val = np.nanmean(vals)
            else:
                mean_val = -np.inf

            if mean_val > best_mean:
                best_mean  = mean_val
                best_probe = probe

        values = []
        for s in sample_ids:
            if s in gene_probe_vals[gene][best_probe]:
                values.append(gene_probe_vals[gene][best_probe][s])
            else:
                values.append(np.nan)

        gene_expr[gene] = pd.Series(values, index=sample_ids)

    expr_wide            = pd.DataFrame(gene_expr).T
    expr_wide.index.name = "gene"

    found   = []
    missing = []
    for gene in TARGET_GENES:
        if gene in expr_wide.index:
            found.append(gene)
        else:
            missing.append(gene)

    found   = sorted(found)
    missing = sorted(missing)

    print("  Genes found: " + str(len(found)) + " / " + str(len(TARGET_GENES)))
    if len(missing) > 0:
        print("  Not in platform: " + str(missing))
    if len(found) > 0:
        expr_wide = expr_wide.loc[found]

    return expr_wide, meta


def to_long(expr, meta):
    """
    Reshape the wide genes-x-samples matrix into long format.

    Wide format has one row per gene and one column per sample.
    Long format has one row per gene-sample combination, which is what
    the downstream DEG analysis scripts expect.
    """
    print("\nReshaping to long format...")

    long = expr.reset_index().melt(
        id_vars="gene",
        var_name="donor_id",
        value_name="cell_density"
    )

    meta_reset = meta.reset_index().rename(columns={"sample_id": "donor_id"})

    keep = ["donor_id"]
    for col in ["diagnosis", "age_years", "age_days", "sex"]:
        if col in meta_reset.columns:
            keep.append(col)

    long = long.merge(meta_reset[keep], on="donor_id", how="left")

    if "age_days" not in long.columns and "age_years" in long.columns:
        long["age_days"] = (long["age_years"] * 365.25).round()

    long["structure_name"]    = "dorsolateral prefrontal cortex"
    long["structure_acronym"] = "DLPFC"
    long["hemisphere"]        = "unknown"

    print("  Shape: " + str(long.shape))
    return long


def save(long, wide, meta):
    long.to_csv(DATA_DIR / "allen_scz_raw.csv",       index=False)
    wide.to_csv(DATA_DIR / "geo_expression_wide.csv")
    meta.to_csv(DATA_DIR / "geo_metadata.csv")

    print("\n" + "=" * 60)
    print("SUCCESS")
    print("=" * 60)
    print("  data/allen_scz_raw.csv   (" + str(len(long)) + " records)")
    print("  Genes   : " + str(long["gene"].nunique()))
    print("  Donors  : " + str(long["donor_id"].nunique()))

    if "diagnosis" in long.columns:
        deduped = long.drop_duplicates("donor_id")
        counts  = deduped["diagnosis"].value_counts()
        for label in counts.index:
            print("  " + str(label) + ": " + str(counts[label]) + " donors")

    print("\nDataset: " + GEO_ACCESSION + " - DLPFC SCZ vs Control (Lanz et al. 2019)")
    print("\nNext: python scripts/02_preprocess.py")


def main():
    print("=" * 60)
    print(" GEO Fetcher - " + GEO_ACCESSION)
    print(" DLPFC Schizophrenia vs Control")
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