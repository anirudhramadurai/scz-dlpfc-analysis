#!/bin/bash
# run_all.sh
# Runs the full scz-dlpfc-analysis pipeline in order.
# Run this from the root of the repository after activating your virtual environment.
#
# Usage:
#   chmod +x run_all.sh
#   ./run_all.sh
#
# Prerequisites:
#   python -m venv venv
#   source venv/bin/activate
#   pip install -r requirements.txt

set -e  # Stop immediately if any script fails

echo "============================================================"
echo " scz-dlpfc-analysis: Full Pipeline"
echo "============================================================"
echo ""

echo "Step 1: Downloading GSE53987 from NCBI GEO (~75 MB, ~5 min)..."
python scripts/01_fetch_geo.py
echo ""

echo "Step 2: Preprocessing (log2 transform, z-score; batch correction skipped - single cohort in DLPFC-only data)..."
python scripts/02_preprocess.py
echo ""

echo "Step 3: Cluster analysis and co-expression figures..."
python scripts/03_cluster_analysis.py
echo ""

echo "Step 4: Differential expression and pathway enrichment..."
python scripts/04_pathway_analysis.py
echo ""

echo "============================================================"
echo " Pipeline complete."
echo " Figures  -> figures/"
echo " Data     -> data/"
echo " Results  -> output/"
echo "============================================================"
