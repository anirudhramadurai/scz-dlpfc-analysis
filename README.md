# scz-dlpfc-analysis

**Bioinformatics analysis of GABAergic interneuron gene expression in the dorsolateral prefrontal cortex in schizophrenia**

*BIME 534 — Biomedical Informatics | University of Utah*

---

## Research Question

Do the expression profiles of GABAergic interneuron subtype markers (parvalbumin+, somatostatin+, VIP+/calretinin+) in postmortem dorsolateral prefrontal cortex (DLPFC) cluster by cell-type identity, by diagnosis (schizophrenia vs. control), or by both — and does this pattern differ from non-interneuron genes included in the same panel?

The secondary question is whether genes that are differentially expressed between schizophrenia and control donors are enriched in specific biological pathways (e.g., GABAergic synapse, glutamate receptor signalling, dopamine pathway).

---

## Biological Background

Schizophrenia affects approximately 1% of the population and is characterized by positive symptoms (hallucinations, delusions), negative symptoms (apathy, avolition), and cognitive deficits (working memory, executive function). The cognitive and negative symptoms are most resistant to treatment and most strongly tied to DLPFC dysfunction [1].

Among the most replicated molecular findings in postmortem schizophrenia brain tissue is a selective reduction in GABAergic cortical interneuron markers, particularly in the DLPFC [1, 2]. The three major non-overlapping interneuron populations — parvalbumin+ (PV+), somatostatin+ (SST+), and calretinin+ (CR+) — show differential vulnerability:

- **PV+ neurons** (basket and chandelier cells, targeting soma and axon initial segments): most studied; most evidence supports decreased PVALB mRNA and PV protein expression, with some studies reporting reduced cell density [1].
- **SST+ neurons** (Martinotti cells, targeting apical dendrites): decreased SST mRNA expression is one of the most consistently replicated findings in DLPFC in schizophrenia [1, 2].
- **CR+ neurons** (targeting other interneurons): largely unaffected in cerebral cortex; this selective sparing is itself a key feature of the interneuron pathology [1].

Reduced GAD1 (GAD67) mRNA expression — encoding the primary enzyme for GABA synthesis — is also consistently found across multiple cortical regions [1, 2]. Critically, pan-interneuron markers (GAD1, GAD2, SLC6A1) show no decrease in cell density, suggesting the pathology reflects a loss of gene expression in surviving cells rather than frank cell loss [2].

These interneuron changes are consistent with GABA, glutamate (NMDA hypofunction), and dopamine hypotheses of schizophrenia, and align with the neurodevelopmental multiple-hit model: SST+ and PV+ neurons originate from the medial ganglionic eminence (MGE), making them specifically vulnerable to hits during a shared developmental window, whereas CR+ neurons arise from the caudal ganglionic eminence and are largely spared [1].

---

## Data Source

Gene expression data were obtained from NCBI GEO accession **GSE53987** [3].

| Field | Value |
|---|---|
| GEO Accession | GSE53987 |
| Platform | GPL570 — Affymetrix Human Genome U133 Plus 2.0 Array |
| Tissue | Postmortem prefrontal cortex (Brodmann Area 46), striatum, hippocampus |
| Groups | Schizophrenia, bipolar disorder, major depressive disorder, matched controls |
| n per group | 19 |
| Total samples | 205 |
| Submission date | January 10, 2014 |
| Deposited by | Thomas A. Lanz, Pfizer Multi-Omics & Biomarkers |

This analysis focuses on the **prefrontal cortex (Brodmann Area 46) samples from schizophrenia and control donors** (48 schizophrenia and 55 control PFC-labelled samples; multiple samples per donor are present across brain regions in the full dataset). Striatum and hippocampus samples are excluded.

**Note on data source selection:** The originally intended data source was the Allen Brain Atlas Human ISH Schizophrenia Study, which used in situ hybridization (ISH) to measure expression of 58 genes at cellular resolution in DLPFC from 19 schizophrenia donors and 33 controls, as described in Guillozet-Bongaarts et al. (2014) [2] and the Allen Human Brain Atlas [4]. The Allen Brain Atlas ISH REST API (api.brain-map.org) is no longer operational; data have migrated to the Brain Knowledge Platform under a different access model. GSE53987 was selected as the replacement data source because it covers the same brain region, the same diagnostic comparison, and includes an overlapping gene set, with the additional advantage that Affymetrix microarray provides continuous quantitative signal for all probes simultaneously, rather than ISH cell density estimates per gene.

---

## Target Gene Panel

Genes were selected to match the published 58-gene panel from Guillozet-Bongaarts et al. (2014) [2], organised into three categories:

**GABAergic interneuron subtype markers**
`GAD1`, `GAD2`, `SLC6A1`, `CALB1`, `CALB2`, `PVALB`, `SST`, `VIP`, `NPY`, `GRIK1`, `ERBB4`, `TAC1`, `TAC3`

**Schizophrenia candidate / associated genes**
`AKT1`, `ARC`, `BDNF`, `CAMK2A`, `CHRNA7`, `CNR1`, `COMT`, `DISC1`, `DLG4`, `FEZ1`, `GRIK4`, `KCNH2`, `NDEL1`, `PAFAH1B1`, `PRODH`, `RELN`, `RGS4`, `SLC1A2`, `SLC1A3`, `TAC1`

**Laminar / cell-type context markers**
`CARTPT`, `CUX2`, `MBP`, `NEFH`, `NR4A2`, `PCP4`, `SYNPR`, `NTNG2`, `CTGF`, `CLDN5`, `CNP`, `B3GALT2`, `C8orf79`, `MFGE8`, `NDNF`, `NOS1AP`, `RORB`, `SCN4B`, `SNCG`, `VAMP1`, `VATIL`

---

## Analysis Pipeline

| Step | Script | Method | Output |
|---|---|---|---|
| 1. Data retrieval | `01_fetch_geo.py` | Direct HTTPS download of GSE53987 SOFT file with custom parsing (GEOparse skipped due to FTP connection issues) | `data/geo_raw.csv`, `data/metadata.csv` |
| 2. Preprocessing | `02_preprocess.py` | RMA-normalized probe values → log₂ + z-score per gene, probe-to-gene mapping, target gene filtering | `data/expression_matrix.csv` |
| 3. Cluster analysis | `03_cluster_analysis.py` | Hierarchical clustering (Ward's linkage, Euclidean), PCA, biclustered gene heatmap, interneuron co-expression matrices, cluster validation (silhouette, ARI, cophenetic r) | `figures/fig1_donor_dendrogram.png`, `fig2_pca.png`, `fig3_gabaergic_coexpression.png`, `fig4_gene_heatmap.png`, `fig5_cluster_validation.png` |
| 4. Differential expression & pathway analysis | `04_pathway_analysis.py` | Welch's t-test + Benjamini–Hochberg FDR correction; over-representation analysis (Fisher's exact test) against KEGG/GO gene sets | `output/differential_expression.csv`, `figures/fig6_volcano.png`, `fig7_barplot_fc.png`, `fig8_pathway_enrichment.png` |

The three-generation framework of pathway analysis (over-representation analysis → functional class scoring → pathway topology) is described in Khatri et al. (2012) [5]. This pipeline implements ORA as the primary enrichment method, appropriate given the moderate sample size and the goal of interpretable pathway-level summaries. Batch effects — a major concern in multi-site postmortem microarray studies — are assessed using principal component analysis as a diagnostic step prior to clustering, following Leek et al. (2010) [6].

---

## Project Structure

```
scz-dlpfc-analysis/
│
├── README.md
├── requirements.txt
├── .gitignore
│
├── scripts/
│   ├── 01_fetch_geo.py            # Download GSE53987 from NCBI GEO via direct HTTPS + SOFT parsing
│   ├── 02_preprocess.py           # Normalize, map probes, filter target genes
│   ├── 03_cluster_analysis.py     # Hierarchical clustering, PCA, co-expression
│   └── 04_pathway_analysis.py     # Differential expression + pathway enrichment
│
├── data/                          # Generated by scripts (gitignored except .gitkeep)
│   └── .gitkeep
│
├── figures/                       # Generated figures (gitignored except .gitkeep)
│   └── .gitkeep
│
└── output/                        # Summary tables (gitignored except .gitkeep)
    └── .gitkeep
```

---

## Setup & Usage

**Requirements:** Python 3.9+

```bash
# 1. Clone the repository
git clone https://github.com/yourusername/scz-dlpfc-analysis.git
cd scz-dlpfc-analysis

# 2. Create a virtual environment (recommended)
python -m venv venv
source venv/bin/activate          # macOS / Linux
venv\Scripts\activate             # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Download data from NCBI GEO (requires internet, ~5–10 min)
python scripts/01_fetch_geo.py

# 5. Run the full analysis pipeline in order
python scripts/02_preprocess.py
python scripts/03_cluster_analysis.py
python scripts/04_pathway_analysis.py
```

All figures are written to `figures/`. All processed data files are written to `data/`. Summary tables are written to `output/`.

---

## References

1. Prkačin MV, Banovac I, Petanjek Z, Hladnik A. Cortical interneurons in schizophrenia – cause or effect? *Croatian Medical Journal*. 2023;64:110–122. doi:10.3325/cmj.2023.64.110

2. Guillozet-Bongaarts AL, Hyde TM, Dalley RA, Hawrylycz MJ, Henry A, Hof PR, Hohmann J, Jones AR, Kuan CL, Royall J, Shen E, Swanson B, Zeng H, Kleinman JE. Altered gene expression in the dorsolateral prefrontal cortex of individuals with schizophrenia. *Molecular Psychiatry*. 2014;19:478–485. doi:10.1038/mp.2013.30. PMID: 23528911

3. Lanz TA, Reinhart V, Sheehan MJ, Rizzo SJS, et al. Postmortem transcriptional profiling reveals widespread increase in inflammation in schizophrenia: a comparison of prefrontal cortex, striatum, and hippocampus among matched tetrads of controls with subjects diagnosed with schizophrenia, bipolar or major depressive disorder. *Translational Psychiatry*. 2019;9(1):151. doi:10.1038/s41398-019-0494-6. PMID: 31123247. GEO: GSE53987.

4. Hawrylycz MJ, Lein ES, Guillozet-Bongaarts AL, Shen EH, Ng L, Miller JA, et al. An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature*. 2012;489:391–399. doi:10.1038/nature11405

5. Khatri P, Sirota M, Butte AJ. Ten years of pathway analysis: current approaches and outstanding challenges. *PLoS Computational Biology*. 2012;8(2):e1002375. doi:10.1371/journal.pcbi.1002375

6. Leek JT, Scharpf RB, Corrada Bravo H, Simcha D, Langmead B, Johnson WE, Geman D, Baggerly K, Irizarry RA. Tackling the widespread and critical impact of batch effects in high-throughput data. *Nature Reviews Genetics*. 2010;11:733–739. doi:10.1038/nrg2825

---

## Acknowledgements

Course: BIME 534, Biomedical Informatics, University of Utah. Data obtained from NCBI Gene Expression Omnibus (GEO), a public repository maintained by the National Center for Biotechnology Information (NCBI). Allen Human Brain Atlas data accessible at human.brain-map.org, produced by the Allen Institute for Brain Science.
