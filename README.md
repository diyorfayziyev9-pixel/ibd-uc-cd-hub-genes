# ibd-uc-cd-hub-genes

R scripts for an in silico differential expression analysis of ulcerative colitis (UC) and Crohn's disease (CD).

## Overview
- Discovery cohort: GSE75214 (UC colon biopsies vs colon controls; CD ileal biopsies vs ileal controls)
- Differential expression: limma (BH-adjusted P < 0.05, |log2FC| > 1)
- Overlap of UC and CD DEGs, GO Biological Process and KEGG enrichment (clusterProfiler)
- PPI network (STRING, combined score >= 0.7) and hub genes (igraph degree)
- Validation: GSE206171 (UC), GSE66407 (CD ileum); Wilcoxon test with BH correction, ROC analysis
- Immune context: xCell enrichment scores and Spearman correlations
- Drug-gene links: DGIdb (hypotheses only)

## Usage
Run the scripts in numerical order (01 to 50). Scripts 45-50 contain the robustness analyses (sensitivity of hub selection, alternative centrality measures, resampling stability), the PTGS2 supplementary validation and the colonic Crohn's disease validation (GSE20881). Summary tables are in the `results` folder.

## Software
R with limma, clusterProfiler, igraph, pROC, xCell, ggplot2 and patchwork. See the Methods section of the article for versions.

## License
MIT
