# ibd-uc-cd-hub-genes

R scripts for an in silico differential expression analysis of ulcerative colitis (UC) and Crohn's disease (CD).

## Overview
- Discovery cohort: GSE75214 (UC colon biopsies vs colon controls; CD ileal biopsies vs ileal controls)
- Differential expression: limma (BH-adjusted P < 0.05, |log2FC| > 1)
- Overlap of UC and CD DEGs, GO Biological Process and KEGG enrichment (clusterProfiler)
- PPI network (STRING, combined score >= 0.7) and hub genes (igraph; genes in the top 10 by degree, betweenness and closeness)
- Robustness of hub selection: nine combinations of log2FC threshold and STRING score, six centrality measures, 200 stratified subsamples
- Validation: GSE206171 (UC), GSE66407 (CD ileum), GSE20881 (colonic CD); Wilcoxon test with BH correction, ROC analysis
- Immune context: xCell enrichment scores and Spearman correlations (BH-adjusted)
- Drug-gene links: DGIdb (hypotheses only)

## Data
All expression data are public in the NCBI Gene Expression Omnibus: GSE75214, GSE206171, GSE66407 and GSE20881. The scripts download or read them as described in `01_setup.R` to `03_select_samples.R` and `19_validation_load.R`, `21_extract_cel.R`, `30_cd_validation_load.R`.

## Usage
Run the scripts in numerical order (01 to 55).
- Scripts 45-50: robustness analyses (sensitivity of hub selection, resampling stability, alternative centrality measures), supplementary validation of PTGS2, and validation in colonic CD (GSE20881).
- Scripts 44, 51, 52 and 53 assemble the figures. Run 44 before 52, because 52 overwrites Figure 7 and S1 with the final versions. Figure 1 (workflow) is drawn by `53_fig1_workflow.py` (Python, matplotlib).
- `54_reference_lookup.R` checks reference metadata in Crossref and the PubMed links of the GEO series. It is not part of the analysis and is only needed to verify the reference list.
- `33_specific_enrichment.R` produces the enrichment plots of UC-specific and CD-specific DEGs (Figure 8).
- `55_figure_fixes.R` redraws the PCA, volcano and heatmap panels of Figures 2 and 3 with larger text, the same group order in all panels (controls first) and gene-level counts in the volcano plots, then calls `44_make_figures.R` to reassemble the figures. It needs the objects created by the earlier scripts (PCA scores, DEG tables, expression matrices), so run it in the same R session after them.

## Output
- `results/`: summary tables (CSV). The supplementary tables of the article come from `sensitivity_grid.csv` (Table S1), `hub_methods_comparison.csv` (Table S2), `resampling_hub_stability.csv` (Table S3) and `validation_CD_colon_GSE20881*.csv` (Table S4)
- `results/figures/`: final figures (Figure 1-8 and Figure S1-S3)

## Software
R 4.6.1 with limma, clusterProfiler, igraph, pROC, GEOquery, xCell, ggplot2, patchwork and pheatmap; Python with matplotlib for Figure 1. See the Methods section of the article for versions.

## Status
This repository accompanies a manuscript in preparation. All results are based on public data and are hypotheses that need experimental confirmation.

## License
MIT
