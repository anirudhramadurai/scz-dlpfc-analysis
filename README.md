# scz-dlpfc-analysis

**Bioinformatics analysis of GABAergic interneuron gene expression in the dorsolateral prefrontal cortex in schizophrenia**

---

## What is this project?

This is an independent computational genomics project asking a specific question about schizophrenia at the molecular level: do the different subtypes of inhibitory neurons in the prefrontal cortex break down together or independently in the disease?

To answer it, I downloaded publicly available postmortem brain gene expression data, built a four-script Python analysis pipeline from scratch, and applied unsupervised machine learning, differential expression analysis, and pathway enrichment to real patient data.

The full pipeline runs in about 15 minutes on any laptop, from raw data download to final figures.

---

## The Biology (plain language first)

### What is schizophrenia, and why the prefrontal cortex?

Schizophrenia is a serious psychiatric disorder affecting about 1% of the global population. It is characterised by three categories of symptoms:

- **Positive symptoms** (things that are added): hallucinations, delusions, disorganised thinking
- **Negative symptoms** (things that are reduced): blunted affect, loss of motivation, social withdrawal
- **Cognitive symptoms**: impaired working memory, difficulty planning, reduced executive function

The cognitive and negative symptoms are the hardest to treat and are most strongly linked to dysfunction of the **dorsolateral prefrontal cortex (DLPFC)**, a brain region located behind the forehead that is critical for reasoning, working memory, and decision-making.

### What are GABAergic interneurons, and why do they matter?

The brain has two broad categories of neurons: excitatory neurons (which activate other neurons) and inhibitory neurons (which suppress activity). Most inhibitory neurons use **GABA** (gamma-aminobutyric acid) as their neurotransmitter. These are called GABAergic neurons.

A specific class of GABAergic neurons called **interneurons** acts as the brain's local circuit regulators. They control the timing and synchrony of excitatory neuron firing, which is essential for the coordinated neural activity that underlies cognition. Without properly functioning interneurons, excitatory circuits become dysregulated, which is thought to contribute to the cognitive symptoms of schizophrenia.

There are several distinct subtypes of GABAergic interneurons, each with a different shape, location, and role in the circuit:

- **PV+ neurons (parvalbumin-positive)**: Fast-spiking basket and chandelier cells. They target the cell body and axon initial segment of excitatory neurons, giving them powerful control over whether those neurons fire at all. Basket cells wrap around the cell body; chandelier cells target the axon initial segment, the very site where action potentials are initiated.
- **SST+ neurons (somatostatin-positive)**: Martinotti cells that extend long processes up to the top layers of the cortex, targeting the apical dendrites of excitatory neurons. They regulate inputs arriving from other brain regions.
- **VIP+ neurons (VIP-positive)**: A third class that primarily targets *other interneurons* rather than excitatory cells directly. By inhibiting inhibitory neurons, VIP+ cells can effectively disinhibit excitatory circuits.
- **CR+ neurons (calretinin-positive)**: Also largely target other interneurons. Generally thought to be spared in schizophrenia, which is itself a key feature of the pathology.

### What does the postmortem literature say?

Decades of postmortem studies consistently find reduced expression of GABAergic interneuron markers in the DLPFC of people with schizophrenia. Key findings include:

- Reduced **GAD1** (GAD67) mRNA, encoding the primary enzyme that synthesises GABA from glutamate. This is one of the most replicated molecular findings in schizophrenia research.
- Reduced **PVALB** (parvalbumin) expression, particularly in deep cortical layers.
- Reduced **SST** (somatostatin) mRNA, also highly replicated.
- Importantly, pan-interneuron markers such as GAD1 and GAD2 show no reduction in *cell density*, which means the cells are still there. The problem appears to be a downregulation of gene expression within surviving cells, not cell death [2].

### The open question this project addresses

If PV+, SST+, and VIP+ interneurons are all affected in schizophrenia, does this happen in a coordinated way (all subtypes going down together as part of a shared pathological process), or do the subtypes show independent patterns of dysregulation (each going down for its own reasons, possibly at different stages or through different mechanisms)?

This distinction matters for understanding whether schizophrenia involves a general failure of GABAergic inhibition, or something more specific to particular circuit elements.

---

## Research Questions

**Primary:** Do the expression profiles of GABAergic interneuron subtype markers (PV+, SST+, VIP+) in postmortem DLPFC co-vary as a coordinated block in schizophrenia, or do subtypes show independent patterns of dysregulation?

**Secondary:** Does unsupervised hierarchical clustering of the gene panel recover donor diagnosis (schizophrenia vs. control) without using the diagnostic label as input, and are differentially expressed genes enriched in specific interneuron subtype categories?

---

## Data

Gene expression data were obtained from **NCBI GEO accession GSE53987** [3].

GEO (Gene Expression Omnibus) is a public repository maintained by NCBI where researchers deposit their gene expression datasets. Anyone can download and reanalyse the data for free, which is the foundation of reproducible science in genomics.

| Field | Value |
|---|---|
| GEO Accession | GSE53987 |
| Platform | GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array) |
| Tissue | Postmortem prefrontal cortex (Brodmann Area 46), striatum, hippocampus |
| Groups | Schizophrenia, bipolar disorder, major depressive disorder, matched controls |
| Total samples | 205 |
| Used in this analysis | 48 SCZ + 55 Control PFC samples (n = 103) |

The dataset uses **Affymetrix microarrays**, a technology that measures the expression level of thousands of genes simultaneously by hybridising labelled RNA to short DNA probes on a chip. Each probe has a known sequence that binds its complementary target RNA; the fluorescence intensity indicates how much of that RNA is present in the sample.

**Note on data source:** This project was originally designed around the Allen Brain Atlas Human ISH Schizophrenia Study (Guillozet-Bongaarts et al. 2014 [2]), which profiled 58 genes by in situ hybridization in DLPFC tissue from 19 SCZ donors and 33 controls. In situ hybridization (ISH) is a technique that labels RNA directly in tissue sections, allowing you to see exactly which cells express a gene and at what density. The Allen Brain Atlas ISH REST API is no longer operational. GSE53987 was selected as the replacement because it covers the same brain region, the same diagnostic comparison, and an overlapping gene set, and Affymetrix microarray provides continuous quantitative signal across all probes simultaneously.

**Gene panel:** 48 of the 58 target genes from Guillozet-Bongaarts et al. were present on the GPL570 platform and used in analysis. Genes absent from the platform (including CHRNA7 and PRODH) were excluded.

---

## Pipeline Overview

Four Python scripts run in sequence. Every output is fully reproducible from the raw GEO download.

```
01_fetch_geo.py   -->   02_preprocess.py   -->   03_cluster_analysis.py   -->   04_pathway_analysis.py
  Download               Filter, batch-           Clustering, PCA,               Differential
  GSE53987               correct, normalize        co-expression                  expression,
  from NCBI GEO          expression data           heatmaps                       enrichment
```

### What each step does

**Step 1: Data retrieval**
Downloads the full GSE53987 SOFT file directly from NCBI via HTTPS. The SOFT format is a plain-text file containing all sample metadata and expression values for the dataset. A custom parser extracts expression values and maps Affymetrix probe IDs to HGNC gene symbols using the platform annotation table. For genes with multiple probes, the probe with the highest mean expression across all samples is selected.

**Step 2: Preprocessing**
Filters to schizophrenia and control donors only (excluding bipolar disorder and MDD, which are present in the full dataset). Applies linear batch correction, log2 transforms the expression values, and z-score normalises per gene across donors.

*Batch correction:* Before correcting, PCA of the raw data showed that the first principal component (explaining 42% of variance) corresponded to processing cohort (samples beginning GSM1304xxx vs GSM1305xxx) rather than diagnosis. This is a technical artefact from samples being processed at different times, a common problem in multi-cohort postmortem datasets. Linear regression was used to remove this systematic shift between cohorts while preserving within-cohort biological variation, following Leek et al. (2010) [5].

*Log2 transformation:* Gene expression values are log2-transformed before analysis. Raw microarray intensity values are right-skewed (a few very highly expressed genes pull the distribution). Log2 transformation compresses the range, making the data more normally distributed and making fold-changes interpretable: a difference of 1 unit on the log2 scale corresponds to a twofold change in expression.

*Z-scoring:* Each gene is normalised to have mean = 0 and standard deviation = 1 across all donors. This puts all genes on a common scale so that clustering is not dominated by highly expressed genes.

**Step 3: Cluster analysis**
Applies Ward's linkage hierarchical clustering to donors, PCA, GABAergic co-expression matrices, and cluster validation metrics (silhouette score, Adjusted Rand Index, cophenetic correlation coefficient).

**Step 4: Differential expression and pathway enrichment**
Welch's t-test per gene with Benjamini-Hochberg FDR correction, followed by over-representation analysis using Fisher's exact test against manually curated interneuron subtype gene sets.

---

## Results

### Does unsupervised clustering recover diagnosis?

**Short answer: No. This is the expected and honest result.**

![Donor hierarchical clustering dendrogram](figures/fig1_donor_dendrogram.png)

This dendrogram shows the result of hierarchical clustering of all 103 donors based solely on their gene expression profiles, with no knowledge of diagnosis. Each leaf at the bottom is one donor; the colour bar shows their actual diagnosis (red = schizophrenia, blue = control). The height at which two branches merge indicates how different those donors are from each other.

If gene expression reliably distinguished schizophrenia from controls, we would expect to see two large clusters, one mostly red and one mostly blue. Instead, the colours are intermixed throughout. The Adjusted Rand Index (ARI = -0.008) quantifies this: ARI of 1.0 would mean perfect recovery of diagnosis; ARI near 0 means the clustering assignments are no better than random with respect to diagnosis.

This is consistent with a well-established finding in postmortem brain transcriptomics: technical variables (RNA quality, postmortem interval) and biological covariates (age, medication history, smoking) typically explain more expression variance than disease status alone. Reporting this honestly is important. It is not a pipeline failure; it is an accurate result about the nature of the data.

![PCA of donor expression profiles](figures/fig2_pca.png)

Principal Component Analysis (PCA) is a dimensionality reduction technique that finds the directions of greatest variance in the data. Here, each dot is a donor projected onto the first two principal components. PC1 captures 42% of total variance. If diagnosis were the dominant signal, SCZ (red) and Control (blue) donors would separate along one of these axes. They do not: the two groups are completely interleaved, confirming the clustering result.

![Cluster validation](figures/fig5_cluster_validation.png)

Cluster validation across k = 2 to 5 clusters. The left panel shows silhouette score (internal validity: how well-separated the clusters are from each other, regardless of diagnosis). The right panel shows ARI vs diagnosis (external validity: how well the clusters align with known diagnosis labels). Internal validity is modest but real. External validity is near zero across all values of k, confirming that the cluster structure in this data does not correspond to diagnostic groupings.

---

### What happens to GABAergic co-expression structure?

**This is the most interesting result in the project.**

![GABAergic interneuron co-expression: SCZ vs Control](figures/fig3_gabaergic_coexpression.png)

This figure shows Pearson correlation matrices for all pairwise combinations of GABAergic interneuron markers, computed separately for schizophrenia donors (left) and control donors (right). Each cell shows the correlation coefficient between two genes across all donors in that group. Red = positive correlation (the two genes tend to go up and down together across donors); blue = negative correlation (when one is high, the other tends to be low).

Gene labels are coloured by subtype: purple = PV+, blue = SST+, green = VIP+, yellow = CB+, orange = Pan-GABA.

**In controls (right panel):** PVALB (purple, PV+) is relatively independent of the other markers. Its correlations with CALB1, NPY, and PENK range from roughly -0.30 to 0.35, close to zero. This makes biological sense: PV+ and SST+ interneurons serve distinct functional roles, and there is no strong reason for their expression levels to be tightly coupled across healthy individuals.

**In schizophrenia (left panel):** that independence is lost. PVALB now correlates strongly with VIP (r = 0.86), GAD1 (r = 0.78), GAD2 (r = 0.83), and SST (r = 0.71). All markers co-vary as a single module.

**What this means:** In healthy cortex, PV+ and SST+ interneurons maintain relatively independent expression profiles, reflecting their distinct roles in regulating cortical circuits. In schizophrenia, this functional separation collapses: all GABAergic markers rise and fall together across donors. This pattern is consistent with coordinated rather than subtype-specific dysregulation, and it directly addresses Research Question 1.

---

### Which specific genes are differentially expressed?

![Gene expression heatmap](figures/fig4_gene_heatmap.png)

This heatmap shows batch-corrected, z-scored expression for all 48 genes across all 103 donors. Donors are sorted by diagnosis (Control on the left, Schizophrenia on the right); the dashed vertical line marks the boundary. Genes are clustered by Ward's linkage along the y-axis. Gene labels on the left are coloured by subtype. The colour of each cell represents the z-score for that donor-gene combination: red = higher than average expression for that gene, blue = lower than average.

A subtle but consistent shift toward blue (lower expression) in the SCZ half is visible for several GABAergic markers, particularly NPY, SST, PENK, PVALB, and GAD1.

![Volcano plot](figures/fig6_volcano.png)

A volcano plot is a standard way to visualise differential expression results. The x-axis shows the log2 fold-change (how much expression differs between schizophrenia and control donors: negative = lower in SCZ, positive = higher in SCZ). The y-axis shows the -log10(p-value): higher points have smaller p-values and are more statistically significant. The dashed horizontal line marks the FDR threshold; points above it passed correction.

Coloured points are significant. Grey points did not pass the thresholds. Gene labels are coloured by subtype.

15 genes were significant at FDR < 0.1, |log2FC| >= 0.1. Every significant GABAergic marker is on the left side (downregulated in schizophrenia). CLDN5, a gene encoding claudin-5 (a tight junction protein in blood-brain barrier endothelial cells), is the only significantly upregulated gene, a finding reported in other schizophrenia postmortem studies.

*Note on thresholds:* The thresholds used here (FDR < 0.1, |log2FC| >= 0.1) are intentionally more permissive than typical genome-wide standards (FDR < 0.05, |log2FC| >= 0.5 or 1.0). This is justified because we are testing only 48 curated genes rather than ~20,000, so the multiple testing burden is far lower. All threshold choices are disclosed.

![Fold-change bar chart](figures/fig7_barplot_fc.png)

The top 20 genes ranked by absolute effect size (|log2FC|), sorted from most negative at the bottom to most positive at the top. Stars (filled diamonds on left edge) mark genes significant after FDR correction. Raw p-values are shown in the right margin. Blue bars = lower in schizophrenia; red = higher in schizophrenia. PVALB, NPY, SST, PENK, and GAD1 all show consistent downregulation with good statistical support.

---

### Are specific subtypes enriched among the differentially expressed genes?

![Pathway enrichment by interneuron subtype](figures/fig8_pathway_enrichment.png)

Over-representation analysis (ORA) tests whether a particular category of genes is more represented among the differentially expressed genes than you would expect by chance, given the background DEG rate across the full panel.

The test used here is Fisher's exact test, a standard nonparametric test for contingency tables. For each subtype, we ask: of all the genes in that subtype category, what fraction are DEGs? And is that fraction significantly higher than the background rate?

The dashed vertical line shows the background DEG rate across the full 48-gene panel (31%).

**SST+: 3 out of 3 genes significant, Fisher's exact p = 0.026.** This is the only statistically significant enrichment. All three SST+ markers in the panel (SST, NPY, PENK) were differentially expressed.

**Excitatory neurons: 0 out of 4 genes.** No enrichment, confirming the effect is specific to inhibitory interneurons rather than a general transcriptional change.

This directly answers Research Question 2: differentially expressed genes are specifically enriched in the SST+ interneuron subtype category.

Note: enrichment was tested against the curated subtype labels defined in this pipeline, not against external pathway databases such as KEGG or GO. External database enrichment was not implemented.

---

## Summary of Findings

| Question | Finding |
|---|---|
| Does clustering recover diagnosis? | No (ARI close to 0). Expected result; technical and biological covariates dominate variance. |
| Do GABAergic subtypes co-vary independently? | No. In SCZ, all markers co-vary as a single module. PV+ independence from SST+ is lost in controls. |
| Which genes are differentially expressed? | 15 genes (FDR < 0.1). All significant GABAergic markers are downregulated. CLDN5 is upregulated. |
| Which subtypes are enriched in DEGs? | SST+ (3/3 genes, Fisher p = 0.026). Excitatory neurons are unaffected (0/4 genes). |

---

## Setup and Usage

**Requirements:** Python 3.9+

```bash
# 1. Clone the repository
git clone https://github.com/anirudhramadurai/scz-dlpfc-analysis.git
cd scz-dlpfc-analysis

# 2. Create a virtual environment
python -m venv venv
source venv/bin/activate          # macOS / Linux
venv\Scripts\activate             # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Download data from NCBI GEO (75 MB, approximately 5 minutes)
python scripts/01_fetch_geo.py

# 5. Run the analysis pipeline in order
python scripts/02_preprocess.py
python scripts/03_cluster_analysis.py
python scripts/04_pathway_analysis.py
```

Figures are written to `figures/`. Data files are written to `data/`. Result tables are written to `output/`.

---

## Limitations

- Sample size (n = 103) limits statistical power, particularly after FDR correction over a small gene panel
- The 48-gene panel is curated and not unbiased; enrichment results reflect prior biological knowledge built into the gene selection
- 10 of the 58 target genes from Guillozet-Bongaarts et al. were absent from the GPL570 platform and could not be analysed
- Batch correction was inferred from GEO sample ID prefix; RNA integrity number (RIN), postmortem interval (PMI), and medication exposure were not available as explicit covariates in this dataset
- 75% of SCZ donors in the original Allen study had antipsychotic evidence at time of death; the degree to which medication confounds expression differences cannot be determined
- ORA was conducted against manually curated subtype labels only, not against external pathway databases (KEGG, GO, Reactome)
- Results are correlational; no causal inference is possible from observational postmortem data

---

## References

1. Prkacin MV, Banovac I, Petanjek Z, Hladnik A. Cortical interneurons in schizophrenia - cause or effect? *Croatian Medical Journal*. 2023;64:110-122. doi:10.3325/cmj.2023.64.110

2. Guillozet-Bongaarts AL, Hyde TM, Dalley RA, Hawrylycz MJ, Henry A, Hof PR, Hohmann J, Jones AR, Kuan CL, Royall J, Shen E, Swanson B, Zeng H, Kleinman JE. Altered gene expression in the dorsolateral prefrontal cortex of individuals with schizophrenia. *Molecular Psychiatry*. 2014;19:478-485. doi:10.1038/mp.2013.30. PMID: 23528911

3. Lanz TA, Reinhart V, Sheehan MJ, Rizzo SJS, et al. Postmortem transcriptional profiling reveals widespread increase in inflammation in schizophrenia: a comparison of prefrontal cortex, striatum, and hippocampus among matched tetrads of controls with subjects diagnosed with schizophrenia, bipolar or major depressive disorder. *Translational Psychiatry*. 2019;9(1):151. doi:10.1038/s41398-019-0494-6. PMID: 31123247. GEO: GSE53987.

4. Hawrylycz MJ, Lein ES, Guillozet-Bongaarts AL, et al. An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature*. 2012;489:391-399. doi:10.1038/nature11405

5. Leek JT, Scharpf RB, Corrada Bravo H, et al. Tackling the widespread and critical impact of batch effects in high-throughput data. *Nature Reviews Genetics*. 2010;11:733-739. doi:10.1038/nrg2825

6. Khatri P, Sirota M, Butte AJ. Ten years of pathway analysis: current approaches and outstanding challenges. *PLoS Computational Biology*. 2012;8(2):e1002375. doi:10.1371/journal.pcbi.1002375

---

## Acknowledgements

Developed as an independent project, drawing on methods from BIME 534 (Biology & Informatics, University of Washington). Data from NCBI Gene Expression Omnibus. Allen Human Brain Atlas at human.brain-map.org, produced by the Allen Institute for Brain Science.
