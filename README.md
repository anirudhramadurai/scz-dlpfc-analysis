# scz-dlpfc-analysis

**Gene expression of GABAergic interneurons in the dorsolateral prefrontal cortex in schizophrenia**

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Data](https://img.shields.io/badge/Data-GEO%3AGSE53987-orange)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53987)

| | |
|---|---|
| **Dataset** | GSE53987 (NCBI GEO): 15 schizophrenia vs. 19 control postmortem DLPFC donors (n = 34) |
| **DEGs** | 3 genes at FDR < 0.1: PVALB (down, padj = 0.011), CARTPT (up, padj = 0.030), RGS4 (down, padj = 0.099) |
| **Clustering** | ARI = 0.198 — weak, not diagnostically useful. Diagnosis is not recovered, as expected in postmortem tissue. |
| **Main observation** | The independence between PV+ and SST+ markers seen in controls is reduced in schizophrenia (exploratory; no pair survives multiple-testing correction at n = 34) |

---

## Overview

I used a public postmortem microarray dataset (GSE53987; 15 schizophrenia and 19 control DLPFC samples, n = 34) to ask a fairly narrow question: when GABAergic interneuron subtypes are disrupted in schizophrenia, do they change together or separately? The analysis is a four-script Python pipeline that runs end to end from the raw GEO download.

Unsupervised clustering did not separate patients from controls (ARI = 0.198). That is the usual outcome for postmortem brain data, where tissue-quality and demographic variation account for more of the signal than diagnosis does, and I treat it as a sanity check rather than a disappointment.

The co-expression analysis was more informative. In controls, the PV+ marker PVALB is largely uncorrelated with the SST+ markers, which fits their separate roles in the circuit. In the schizophrenia group that separation is weaker: PVALB correlates more strongly with GAD1, NPY, and SST, and the markers drift toward moving as one block. I want to be careful here — formal Fisher z-tests on all 36 pairwise correlations leave nothing significant after Benjamini–Hochberg correction at n = 34, so this is a hypothesis worth following up, not a settled result.

Three genes are differentially expressed at FDR < 0.1: PVALB (PV+, down, padj = 0.011), CARTPT (up, padj = 0.030), and RGS4 (down, padj = 0.099). No interneuron subtype is enriched among the DEGs, which is unsurprising at this sample size. The two obvious caveats: the panel is 48 curated genes rather than the whole transcriptome, and although RIN, PMI, pH, age, and sex are in the metadata and balanced across groups (all Welch p > 0.35), I did not fold them into the differential-expression model.

---

## What this project is

A reproducible Python pipeline applied to postmortem brain expression data, aimed at one question about schizophrenia: do the inhibitory-neuron subtypes of the prefrontal cortex break down together or independently? The whole thing — download to figures — takes about 15 minutes on a laptop.

---

## Background

### Schizophrenia and the prefrontal cortex

Schizophrenia affects roughly 1% of people worldwide [1]. Its symptoms fall into three groups: positive symptoms such as hallucinations and delusions; negative symptoms such as blunted affect and social withdrawal; and cognitive symptoms such as impaired working memory and planning. The cognitive and negative symptoms are the least treatable, and both point toward the **dorsolateral prefrontal cortex (DLPFC)**, the frontal-lobe region behind reasoning, working memory, and executive control [1].

### GABAergic interneurons

Neurons come in two broad flavors: excitatory (they drive other neurons) and inhibitory (they hold activity down). Most inhibitory neurons signal with **GABA**, so they are called GABAergic. A subset of these, the **interneurons**, act as local traffic controllers — they set the timing and synchrony of excitatory firing, which the brain needs for coherent activity during thought [1, 10]. When interneurons work poorly, excitatory circuits lose their rhythm, and this is one of the leading circuit-level explanations for the cognitive symptoms of schizophrenia [1, 10].

Interneurons are not a single population. The subtypes differ in shape, wiring, and job [1]:

- **PV+ (parvalbumin)** — fast-spiking basket and chandelier cells. Basket cells clasp the cell body of excitatory neurons; chandelier cells sit on the axon initial segment, where action potentials start. Between them they hold strong veto power over whether an excitatory neuron fires.
- **SST+ (somatostatin)** — Martinotti cells that reach up to the superficial cortical layers and target the distal dendrites of excitatory neurons, shaping the long-range input those cells receive.
- **VIP+ (VIP)** — cells that mostly inhibit other interneurons rather than excitatory neurons, so they can release excitatory circuits from inhibition indirectly.
- **CR+ (calretinin)** — also largely target other interneurons and are generally reported as spared in schizophrenia, which itself helps characterize the pathology [1].
- **CB+ (calbindin)** — a marker carried by a subset of interneurons, including some SST+ cells in the DLPFC; used here as a co-expression marker, with substantial overlap with SST+ in layers II–III [1, 10].

PV+ and SST+ interneurons share an origin in the medial ganglionic eminence (MGE), which is one reason a shared developmental insult could affect both; CR+ cells come from the caudal ganglionic eminence (CGE) on a separate trajectory [1].

### What the postmortem literature reports

Postmortem studies have converged, over decades, on reduced expression of GABAergic markers in the schizophrenia DLPFC [1, 2, 10]:

- Reduced **GAD1** (GAD67) mRNA — the enzyme that makes GABA from glutamate, and one of the most reproducible molecular findings in the field [1, 2].
- Reduced **PVALB**, especially in the deeper layers [1, 2].
- Reduced **SST** mRNA, replicated across independent cohorts [1, 2].
- Notably, interneuron *density* (by pan-interneuron markers such as GAD1 at the cellular level) is not reduced [2, 10]. The signal looks like downregulation inside surviving cells rather than cell loss.

Together these support the GABAergic hypothesis: weakened cortical inhibition contributes to prefrontal dysfunction and cognitive symptoms [1, 10].

### The open question

If PV+, SST+, and VIP+ interneurons are all affected, is the disruption coordinated across subtypes or specific to each? Whether the DLPFC shows a broad inhibitory failure or a targeted hit to particular circuit elements remains open, and the answer matters for how we think about the disease [1].

---

## Questions

**Primary.** Do the GABAergic subtype markers (PV+, SST+, VIP+) in the postmortem DLPFC move as one coordinated block in schizophrenia, or do the subtypes vary independently?

**Secondary.** Does unsupervised hierarchical clustering recover diagnosis without being told the labels, and are the differentially expressed genes concentrated in any one interneuron subtype?

---

## Data

Expression data come from **NCBI GEO accession GSE53987** [3]. GEO is NCBI's public repository for expression datasets, where deposited data can be downloaded and reanalyzed freely [8].

| Field | Value |
|---|---|
| GEO accession | GSE53987 |
| Platform | GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array) |
| Tissue | Postmortem prefrontal cortex (Brodmann Area 46), striatum, hippocampus |
| Groups | Schizophrenia, bipolar disorder, major depressive disorder, matched controls |
| Total samples | 205 |
| Used here | 15 SCZ + 19 control DLPFC samples (n = 34) |

The platform is an **Affymetrix microarray**, which reads out thousands of genes at once: labeled RNA hybridizes to short DNA probes of known sequence, and each probe's fluorescence reflects how much of its target RNA was present.

**On the data source.** The project was first built around the Allen Brain Atlas Human ISH Schizophrenia Study [4] (Guillozet-Bongaarts et al. 2014 [2]), which profiled 58 genes by in situ hybridization in DLPFC tissue from 19 SCZ donors and 33 controls. In situ hybridization labels RNA inside intact tissue sections, so you can see which cells express a gene and how densely. That resource's REST API is no longer operational, so I moved to GSE53987, which covers the same region, the same comparison, and an overlapping gene set, and gives continuous quantitative signal across all probes at once.

**Gene panel.** Of the 58 target genes from Guillozet-Bongaarts et al., 50 were carried into the target list; CHRNA7 and PRODH are absent from the GPL570 platform, leaving 48 analyzed. (The remaining 8 of the original 58 were not carried forward.)

---

## Pipeline

Four scripts run in order. Every output regenerates from the raw GEO download.

```
01_fetch_geo.py   -->   02_preprocess.py   -->   03_cluster_analysis.py   -->   04_pathway_analysis.py
  Download               Filter, normalize        Clustering, PCA,               Differential
  GSE53987               (DLPFC only)              co-expression                  expression,
  from NCBI GEO                                    heatmaps                       enrichment
```

**Step 1 — retrieval.** Downloads the GSE53987 SOFT file from NCBI over HTTPS. SOFT (Simple Omnibus Format in Text) is the plain-text container GEO uses for metadata and expression values [8]. A parser pulls the values and maps Affymetrix probe IDs to gene symbols via the platform annotation table; where a gene has several probes, the one with the highest mean expression is kept.

**Step 2 — preprocessing.** Keeps schizophrenia and control DLPFC donors only (bipolar, MDD, hippocampus, and striatum samples are dropped), log2-transforms the intensities, and z-scores each gene across donors.

- *Batch correction.* Once the data are restricted to DLPFC, all 34 donors sit in a single processing cohort (GSM1304xxx), and no batch effect is detectable, so none is applied. The batch effect that showed up in the full 205-sample set was regional biology (DLPFC vs. hippocampus vs. striatum), not a processing artifact within the DLPFC.
- *Log2 transform.* Microarray intensities are right-skewed; log2 compresses the range and makes fold-changes readable, since one log2 unit is a twofold change.
- *Z-scoring.* Each gene is centered to mean 0, SD 1, so that clustering is not dominated by the most highly expressed genes.

**Step 3 — clustering.** Ward's-linkage hierarchical clustering on donors, PCA to look for diagnosis separation, separate GABAergic co-expression matrices for the SCZ and control groups, and cluster validation via silhouette score, Adjusted Rand Index, and cophenetic correlation [9].

**Step 4 — differential expression and enrichment.** Welch's t-test per gene with Benjamini–Hochberg FDR correction [7], then over-representation analysis (Fisher's exact test) against manually curated interneuron subtype gene sets [6].

---

## Results

### Does clustering recover diagnosis?

No, and that is what postmortem data usually does.

![Donor hierarchical clustering dendrogram](figures/fig1_donor_dendrogram.png)

Each leaf is one of the 34 DLPFC donors, clustered on expression alone with no access to the diagnosis label; the color bar underneath shows the true diagnosis (red = schizophrenia, blue = control), and branch height reflects dissimilarity. If expression tracked diagnosis, the tree would split into a mostly-red and a mostly-blue clade. It does not — the colors are interleaved. The Adjusted Rand Index puts a number on this: 1.0 is perfect recovery, 0 is chance, and 0.198 sits close to the low end. There is weak positive alignment, so the clusters are not pure noise, but they do not separate the groups.

This matches a well-documented property of postmortem brain transcriptomics: RNA quality, postmortem interval, age, and medication history usually explain more variance than diagnosis [5]. The clustering behaving this way is a point in the pipeline's favor, not against it.

![PCA of donor expression profiles](figures/fig2_pca.png)

PCA projects the donors onto the two directions of greatest variance [9]. If diagnosis dominated, the red and blue points would fall on opposite sides of one axis; instead they overlap, which agrees with the dendrogram. (An earlier version reported 42% variance on PC1 from a dataset that had accidentally kept hippocampus and striatum samples; this PCA uses DLPFC only.)

![Cluster validation](figures/fig5_cluster_validation.png)

Validation across k = 2 to 5. Silhouette score (left) measures how internally distinct the clusters are [9]; ARI vs. diagnosis (right) measures agreement with the true labels. Internal structure is modest but present; external agreement stays weak throughout (ARI 0.09–0.20).

### What happens to GABAergic co-expression?

![GABAergic interneuron co-expression: SCZ vs Control](figures/fig3_gabaergic_coexpression.png)

These are Pearson correlation matrices over all pairs of GABAergic markers, computed separately for the schizophrenia (left) and control (right) donors. Warm colors are positive correlations, cool colors negative; both panels share the same gene ordering, so a given cell is comparable across the two. Marker labels are colored by subtype (purple = PV+, blue = SST+, green = VIP+, yellow = CB+, orange = pan-GABA).

In controls, PVALB (the PV+ marker) is close to independent of the others — its correlations with CALB1, NPY, and PENK run from about −0.30 to 0.35, i.e. near zero. That independence is what you would expect from cells with distinct jobs [1].

In the schizophrenia group (n = 15), those same PVALB correlations are higher — 0.56 with GAD1, 0.58 with NPY, 0.25 with SST — and the panel as a whole shifts toward positive, with the markers tending to move together [10].

Read one way, this suggests the subtype-specific independence of healthy cortex gives way to a more uniform, block-like pattern in schizophrenia — closer to a broad inhibitory disturbance than a hit confined to one subtype, which is the substance of the primary question. Read carefully, it is provisional: none of the 36 pairwise differences survives multiple-testing correction at this n, so I present it as a lead, not a conclusion.

### Which genes are differentially expressed?

![Gene expression heatmap](figures/fig4_gene_heatmap.png)

Z-scored expression for all 48 genes across all 34 donors, with donors sorted by diagnosis (controls left of the dashed line, schizophrenia right) and genes clustered by Ward's linkage. Red is above the gene's mean, blue below. The clearest shift toward blue on the schizophrenia side is PVALB; SST, NPY, and PENK lean the same way but do not clear FDR at n = 34.

![Volcano plot](figures/fig6_volcano.png)

The x-axis is log2 fold-change (negative = lower in schizophrenia); the y-axis is −log10 of the *adjusted* p-value, so the dashed line at padj = 0.1 is a genuine FDR threshold and points above it are the significant genes. Three clear the line: PVALB (PV+, down) and CARTPT (up) comfortably, and RGS4 (down) right at the boundary (padj = 0.099).

CARTPT is the odd one out. It encodes a neuropeptide (cocaine- and amphetamine-regulated transcript) tied to energy balance and stress response, and its *up*regulation does not fit the general GABAergic-downregulation story — it is worth a closer look rather than an explanation I want to force. The drop from the 15 DEGs an earlier version reported traces to a fixed sample-selection bug: that run had pooled DLPFC, hippocampus, and striatum together.

*On thresholds.* FDR < 0.1 with |log2FC| ≥ 0.1 is looser than the genome-wide norm of FDR < 0.05, |log2FC| ≥ 0.5. That is deliberate: 48 curated genes carry a fraction of the multiple-testing burden of a ~20,000-gene screen, so the looser cut is defensible here [6, 7].

![Fold-change bar chart](figures/fig7_barplot_fc.png)

The 20 genes with the largest absolute log2FC, most-negative at the bottom, with FDR-significant genes starred and raw p-values in the margin. Blue is lower in schizophrenia, red higher. PVALB, NPY, SST, PENK, and GAD1 all point down, in line with Guillozet-Bongaarts et al. (2014) [2].

### Is any subtype enriched among the DEGs?

![Pathway enrichment by interneuron subtype](figures/fig8_pathway_enrichment.png)

Over-representation analysis asks whether a subtype's genes are hit more often than the panel-wide DEG rate would predict [6], using Fisher's exact test on the 2×2 table of (DEG / not) × (in category / not). The dashed line marks the background rate across all 48 genes (6%).

No subtype reaches significance, which follows from the sample size. At the single-gene level, PVALB is the only DEG among the interneuron markers; SST, NPY, and PENK trend down as expected but do not clear FDR. Excitatory neurons contribute 0 of 4 DEGs, consistent with an effect confined to the inhibitory side — the direction the secondary question anticipated, though it is a directional read rather than a statistically supported one at n = 34.

Enrichment here is tested against the curated subtype labels defined in the pipeline, not against external databases such as KEGG or GO [6]; external-database enrichment is not implemented.

---

## Summary

| Question | Finding |
|---|---|
| Does clustering recover diagnosis? | No. ARI = 0.198 — weak positive alignment, not diagnostically useful. Covariates dominate the variance [5]. |
| Do the GABAergic subtypes vary independently? | Less so in schizophrenia: the markers trend toward one block, though no pairwise difference survives correction at n = 34. |
| Which genes are differentially expressed? | Three at FDR < 0.1: PVALB (down), CARTPT (up), RGS4 (down, borderline). Power-limited [7]. |
| Any subtype enriched among DEGs? | None significant. PVALB is the only interneuron-marker DEG [6]. |

---

## Setup

**Requires** Python 3.9+.

```bash
# 1. Clone
git clone https://github.com/anirudhramadurai/scz-dlpfc-analysis.git
cd scz-dlpfc-analysis

# 2. Virtual environment
python -m venv venv
source venv/bin/activate          # macOS / Linux
venv/Scripts/activate             # Windows

# 3. Dependencies
pip install -r requirements.txt

# 4. Run everything
chmod +x run_all.sh
./run_all.sh
```

Or step by step:

```bash
python scripts/01_fetch_geo.py     # download from NCBI GEO (~75 MB, ~5 min)
python scripts/02_preprocess.py
python scripts/03_cluster_analysis.py
python scripts/04_pathway_analysis.py
```

Figures land in `figures/`, data in `data/`, result tables in `output/`.

---

## Limitations

- **Sample size.** With 34 DLPFC donors, power is low; only the largest effects clear FDR [7].
- **Curated panel.** 48 genes chosen from prior biology, not the whole transcriptome — so any enrichment reflects the panel's built-in priors [6]. Of the 58 source genes from Guillozet-Bongaarts et al., 50 entered the target list; CHRNA7 and PRODH are off the GPL570 platform, leaving 48 [2].
- **Covariates not modeled.** No batch correction is needed (single DLPFC cohort). RIN, PMI, pH, age, and sex are in the metadata and balanced across groups (all Welch p > 0.35), but I did not add them to the DE model — putting them into a linear model is the natural next step.

  | Covariate | SCZ (n=15) | Control (n=19) | Welch *p* |
  |---|---|---|---|
  | Age (years) | 46.0 ± 8.6 | 48.1 ± 10.6 | 0.54 |
  | PMI (hours) | 18.9 ± 6.7 | 19.5 ± 5.1 | 0.77 |
  | Brain pH | 6.5 ± 0.4 | 6.6 ± 0.2 | 0.54 |
  | RIN | 7.6 ± 0.7 | 7.8 ± 0.6 | 0.38 |
  | Sex (M/F) | 7/8 | 10/9 | 1.00 (Fisher) |

- **Medication.** Antipsychotic exposure at death is not in the GEO metadata [3], so medication effects cannot be separated from disease effects — a limitation shared across the postmortem field [2].
- **Enrichment scope.** Tested against the curated subtype labels only, not KEGG/GO/Reactome [6].
- **Assay.** Microarray reads relative abundance across many probes at once but is less sensitive than RNA-seq for low-expression genes.
- **Design.** Postmortem data are observational; nothing here supports a causal claim.

## Future directions

- Rerun on RNA-seq cohorts (CommonMind, PsychENCODE) to see whether the pattern holds with a more sensitive assay.
- Add age, PMI, RIN, and medication as covariates where the metadata allow.
- Extend enrichment to GO biological-process terms and single-cell-derived cell-type signatures.
- Compare co-expression structure across the bipolar and MDD groups in GSE53987 to test how specific the schizophrenia pattern is.
- Try WGCNA to find diagnosis-associated modules rather than testing genes one at a time.

---

## FAQ

**If clustering didn't recover diagnosis, did the analysis fail?**
No. Postmortem brain expression is dominated by tissue quality, postmortem interval, age, and medication, so a clean split by diagnosis would be the surprising outcome. The interleaving is what the literature predicts and is a sign the pipeline isn't manufacturing structure that isn't there.

**Why only 3 DEGs when an earlier version had 15?**
The 15-gene version had accidentally pooled DLPFC, hippocampus, and striatum. What looked like a batch effect there was real regional biology. Restricting to DLPFC (n = 34) leaves 3 genes past FDR — fewer, but correct.

**Why FDR < 0.1 rather than 0.05?**
The 0.05 convention is calibrated for ~20,000-gene screens. This panel is 48 genes, a much lighter multiple-testing load, so 0.1 is a reasonable and stated choice.

**What would reduced co-expression independence mean biologically?**
In controls, PV+ and SST+ interneurons do different jobs and their markers vary fairly independently. If that independence weakens and the markers move together, it points toward a broader inhibitory disturbance rather than a subtype-specific one — as a hypothesis, given the sample size.

**PVALB was the only significant interneuron marker — are only PV+ cells affected?**
Not necessarily. PVALB, GAD1, SST, NPY, and PENK all trend down together, matching Guillozet-Bongaarts et al. (2014); the others simply don't clear FDR at n = 34, which is a power issue, not evidence of no effect.

**Could antipsychotics confound the PVALB result?**
Yes — a real limitation. Chronic antipsychotic use changes expression, and without exposure records (absent from GSE53987) it can't be separated out. This applies across the postmortem field.

**Why microarray instead of RNA-seq?**
GSE53987 is the best public postmortem DLPFC set covering this gene panel with balanced groups. RNA-seq cohorts are listed as future work; they detect low-abundance transcripts and splicing differences better.

---

## References

1. Prkacin MV, Banovac I, Petanjek Z, Hladnik A. Cortical interneurons in schizophrenia - cause or effect? *Croatian Medical Journal*. 2023;64:110-122. doi:10.3325/cmj.2023.64.110

2. Guillozet-Bongaarts AL, Hyde TM, Dalley RA, Hawrylycz MJ, Henry A, Hof PR, Hohmann J, Jones AR, Kuan CL, Royall J, Shen E, Swanson B, Zeng H, Kleinman JE. Altered gene expression in the dorsolateral prefrontal cortex of individuals with schizophrenia. *Molecular Psychiatry*. 2014;19:478-485. doi:10.1038/mp.2013.30. PMID: 23528911

3. Lanz TA, Reinhart V, Sheehan MJ, Rizzo SJS, et al. Postmortem transcriptional profiling reveals widespread increase in inflammation in schizophrenia: a comparison of prefrontal cortex, striatum, and hippocampus among matched tetrads of controls with subjects diagnosed with schizophrenia, bipolar or major depressive disorder. *Translational Psychiatry*. 2019;9(1):151. doi:10.1038/s41398-019-0492-8. PMID: 31123247. GEO: GSE53987.

4. Hawrylycz MJ, Lein ES, Guillozet-Bongaarts AL, et al. An anatomically comprehensive atlas of the adult human brain transcriptome. *Nature*. 2012;489:391-399. doi:10.1038/nature11405

5. Leek JT, Scharpf RB, Corrada Bravo H, et al. Tackling the widespread and critical impact of batch effects in high-throughput data. *Nature Reviews Genetics*. 2010;11:733-739. doi:10.1038/nrg2825

6. Khatri P, Sirota M, Butte AJ. Ten years of pathway analysis: current approaches and outstanding challenges. *PLoS Computational Biology*. 2012;8(2):e1002375. doi:10.1371/journal.pcbi.1002375

7. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society Series B*. 1995;57(1):289-300.

8. Barrett T, Wilhite SE, Ledoux P, et al. NCBI GEO: archive for gene expression and epigenomics data sets: 23-year update. *Nucleic Acids Research*. 2024;52(D1):D138-D144. doi:10.1093/nar/gkad965. PMID: 37933855

9. Tan PN, Steinbach M, Kumar V. *Introduction to Data Mining*. Pearson; 2006. Chapter 8: Cluster Analysis.

10. Dienel SJ, Lewis DA. Alterations in cortical interneurons and cognitive function in schizophrenia. *Neurobiology of Disease*. 2019;131:104208. PMC6309598.

---

## Acknowledgements

An independent project, using methods from BIME 534 (Biology & Informatics, University of Washington) — reproducible analysis of public high-throughput data and translational reading of omics results. Data from NCBI Gene Expression Omnibus. Allen Human Brain Atlas at human.brain-map.org, produced by the Allen Institute for Brain Science.