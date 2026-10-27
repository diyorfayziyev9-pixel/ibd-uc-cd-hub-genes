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

## Usage
Run the scripts in numerical order (01 to 53).
- Scripts 45-50: robustness analyses (sensitivity of hub selection, resampling stability, alternative centrality measures), supplementary validation of PTGS2, and validation in colonic CD (GSE20881).
- Scripts 44, 51, 52 and 53 assemble the figures. Run 44 before 52, because 52 overwrites Figure 7 and S1 with the final versions. Figure 1 (workflow) is drawn by `53_fig1_workflow.py` (Python, matplotlib).

## Output
- `results/`: summary tables (CSV), including the supplementary tables (sensitivity grid, centrality comparison, resampling stability, validation tables)
- `results/figures/`: final figures (Figure 1-8 and Figure S1-S3)

## Software
R 4.6.1 with limma, clusterProfiler, igraph, pROC, GEOquery, xCell, ggplot2, patchwork and pheatmap; Python with matplotlib for Figure 1. See the Methods section of the article for versions.

## License
MIT
