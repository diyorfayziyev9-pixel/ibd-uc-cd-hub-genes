library(ggplot2)
library(pheatmap)

# --- Yagona rang qoidasi: guruh nomiga qarab avtomatik ---
group_cols <- function(lv) {
  lv <- as.character(unique(lv))
  cols <- ifelse(grepl("control|normal|healthy|ctrl", lv, ignore.case = TRUE), "#5B7FA6",
                 ifelse(grepl("CD|crohn", lv, ignore.case = FALSE), "#8E3B6E", "#D1495B"))
  setNames(cols, lv)
}

# --- PCA ---
pca_plot <- function(pcdf, title, file) {
  p <- ggplot(pcdf, aes(PC1, PC2, color = group)) +
    geom_point(size = 3, alpha = 0.85) +
    scale_color_manual(values = group_cols(pcdf$group), name = "Group") +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank()) +
    labs(title = title)
  print(p)
  ggsave(file, p, width = 7, height = 6, dpi = 300)
}

# --- Volcano ---
volcano <- function(df, title, file) {
  df$status <- ifelse(df$adj.P.Val < 0.05 & df$logFC > 1, "Up",
                      ifelse(df$adj.P.Val < 0.05 & df$logFC < -1, "Down", "NS"))
  n_up   <- sum(df$status == "Up")
  n_down <- sum(df$status == "Down")
  p <- ggplot(df, aes(logFC, -log10(adj.P.Val), colour = status)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_colour_manual(values = c(Up = "#C0392B", Down = "#2E6F9E", NS = "grey75"),
                        name = NULL) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank()) +
    labs(title = title,
         subtitle = paste0("Up: ", n_up, "   Down: ", n_down),
         x = "log2 fold change", y = "-log10 adj. P")
  print(p)
  ggsave(file, p, width = 7, height = 6, dpi = 300)
}

# --- Heatmap ---
make_heatmap <- function(deg_clean, ex, ids, groups, file, n = 40) {
  top <- head(deg_clean[order(deg_clean$adj.P.Val), ], n)
  m <- ex[as.character(top$ID), ]
  rownames(m) <- top$symbol
  ann <- data.frame(Group = groups, row.names = colnames(m))
  ann_cols <- list(Group = group_cols(ann$Group))
  png(file, width = 2200, height = 2400, res = 300)
  pheatmap(m, scale = "row", annotation_col = ann,
           annotation_colors = ann_cols,
           show_colnames = FALSE, fontsize_row = 7,
           color = colorRampPalette(c("#2E6F9E", "white", "#C0392B"))(100))
  dev.off()
}