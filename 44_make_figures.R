# 44_make_figures.R
# Maqola figuralarini tayyor PNG'lardan yig'ish (rang qoidasi: 43_style_functions.R)
library(ggplot2); library(patchwork); library(png); library(grid)
dir.create("results/figures", showWarnings = FALSE)

make_fig <- function(files, ncol, out, w, h, tag = TRUE) {
  plots <- lapply(files, function(f) {
    img <- readPNG(f)
    a <- nrow(img) / ncol(img)
    ggplot() +
      annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = a, interpolate = TRUE) +
      coord_fixed(xlim = c(0, 1), ylim = c(0, a), expand = FALSE) +
      theme_void() + theme(plot.margin = margin(28, 4, 4, 28))
  })
  fig <- wrap_plots(plots, ncol = ncol) + plot_annotation(tag_levels = if (tag) "A" else NULL) &
    theme(plot.tag = element_text(face = "bold", size = 22))
  ggsave(out, fig, width = w, height = h, dpi = 300)
}

make_fig(c("results/intermediate/PCA_UC_v2.png", "results/intermediate/PCA_CD_v2.png",
           "results/intermediate/volcano_UC_v2.png", "results/intermediate/volcano_CD_v2.png"),
         2, "results/figures/Figure2_PCA_volcano.png", 14, 12)
make_fig(c("results/intermediate/Venn_UC_CD_v2.png", "results/intermediate/heatmap_UC_v2.png", "results/intermediate/heatmap_CD_v2.png"),
         3, "results/figures/Figure3_Venn_heatmaps.png", 16, 6.5)
make_fig(c("results/intermediate/GO_BP_dotplot.png", "results/intermediate/KEGG_dotplot.png"),
         2, "results/figures/Figure4_GO_KEGG.png", 16, 7)
make_fig("results/intermediate/PPI_network_hubs.png",
         1, "results/figures/Figure5_PPI_hubs.png", 9, 7.5, tag = FALSE)
make_fig(c("results/intermediate/validation_boxplot_v2.png", "results/intermediate/validation_ROC_v2.png",
           "results/intermediate/validation_boxplot_CD_v2.png", "results/intermediate/validation_ROC_CD_v2.png"),
         2, "results/figures/Figure6_validation.png", 16, 12)
make_fig(c("results/intermediate/xcell_heatmap.png", "results/intermediate/hub_immune_corr.png"),
         2, "results/figures/Figure7_immune.png", 16, 7)
make_fig(c("results/intermediate/GO_BP_UC_specific.png", "results/intermediate/GO_BP_CD_specific.png"),
         2, "results/figures/Figure8_specific_GO.png", 16, 7)
make_fig(c("results/intermediate/hub_immune_corr_UC.png", "results/intermediate/hub_immune_corr_CD.png"),
         2, "results/figures/FigureS1_within_group_corr.png", 16, 6)
