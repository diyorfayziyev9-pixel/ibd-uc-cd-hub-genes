# 55_figure_fixes.R
# Final fixes for Figures 2 and 3 (readability and consistency):
#   - PCA and volcano plots: larger text (base_size 17), same group order and
#     labels (control, UC / control, CD); volcano counts are gene-level DEGs.
#   - Heatmaps: top 25 genes, larger gene labels, controls always on the left
#     (columns ordered by group, not clustered).
# Requires objects from the earlier scripts: pc, pc_cd (PCA scores with a
# 'group' column), deg_uc_clean, deg_cd_clean, ex_uc, ex_cd, group_uc, group_cd.
# Then Figure 2 and Figure 3 are reassembled by 44_make_figures.R.

library(ggplot2)
library(pheatmap)
source("scripts/43_style_functions.R")   # group_cols()

lab <- function(g, disease) {
  g <- as.character(g)
  factor(ifelse(grepl("control|normal|healthy|ctrl", g, ignore.case = TRUE),
                "control", disease),
         levels = c("control", disease))
}

# --- PCA ---
pca2 <- function(pcdf, disease, file) {
  pcdf$group <- lab(pcdf$group, disease)
  p <- ggplot(pcdf, aes(PC1, PC2, color = group)) +
    geom_point(size = 3, alpha = 0.85) +
    scale_color_manual(values = group_cols(levels(pcdf$group)), name = "Group") +
    theme_bw(base_size = 17) + theme(panel.grid.minor = element_blank()) +
    labs(title = paste0("PCA: ", disease, " vs Control"))
  ggsave(file, p, width = 7, height = 6, dpi = 300)
}
pca2(pc,    "UC", "results/intermediate/PCA_UC_v2.png")
pca2(pc_cd, "CD", "results/intermediate/PCA_CD_v2.png")

# --- Volcano ---
volcano2 <- function(df, disease, file) {
  df$status <- ifelse(df$adj.P.Val < 0.05 & df$logFC > 1, "Up",
                      ifelse(df$adj.P.Val < 0.05 & df$logFC < -1, "Down", "NS"))
  p <- ggplot(df, aes(logFC, -log10(adj.P.Val), colour = status)) +
    geom_point(alpha = 0.6, size = 1.4) +
    scale_colour_manual(values = c(Up = "#C0392B", Down = "#2E6F9E", NS = "grey75"),
                        name = NULL) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    theme_bw(base_size = 17) + theme(panel.grid.minor = element_blank()) +
    labs(title = paste0(disease, " vs Control"),
         subtitle = paste0("Up: ", sum(df$status == "Up"),
                           "   Down: ", sum(df$status == "Down")),
         x = "log2 fold change", y = "-log10 adj. P")
  ggsave(file, p, width = 7, height = 6, dpi = 300)
}
volcano2(deg_uc_clean, "UC", "results/intermediate/volcano_UC_v2.png")
volcano2(deg_cd_clean, "CD", "results/intermediate/volcano_CD_v2.png")

# --- Heatmap (controls on the left, columns not clustered) ---
new_heatmap <- function(deg_clean, ex, groups, file, n = 25) {
  top <- head(deg_clean[order(deg_clean$adj.P.Val), ], n)
  m <- ex[as.character(top$ID), ]
  rownames(m) <- top$symbol
  g <- as.character(groups)
  is_ctrl <- grepl("control|normal|healthy|ctrl", g, ignore.case = TRUE)
  o <- order(!is_ctrl)                      # controls first
  m <- m[, o]; g <- g[o]
  ann <- data.frame(Group = g, row.names = colnames(m))
  png(file, width = 1200, height = 1800, res = 300)
  pheatmap(m, scale = "row", annotation_col = ann,
           annotation_colors = list(Group = group_cols(ann$Group)),
           cluster_cols = FALSE, annotation_names_col = FALSE,
           show_colnames = FALSE, fontsize = 11, fontsize_row = 11,
           treeheight_row = 15, border_color = NA,
           color = colorRampPalette(c("#2E6F9E", "white", "#C0392B"))(100))
  dev.off()
}
new_heatmap(deg_uc_clean, ex_uc, group_uc, "results/intermediate/heatmap_UC_v2.png")
new_heatmap(deg_cd_clean, ex_cd, group_cd, "results/intermediate/heatmap_CD_v2.png")

# --- Reassemble Figures 2 and 3 ---
source("scripts/44_make_figures.R")