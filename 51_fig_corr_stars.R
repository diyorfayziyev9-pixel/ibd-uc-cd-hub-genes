library(ggplot2)

star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
cell_order <- c("Neutrophils","Monocytes","Macrophages M1","Macrophages M2","cDC","aDC",
                "Plasma cells","B-cells","Tregs","CD8+ Tcm","CD4+ memory T-cells")
gene_order <- rev(c("IL6","IL1B","CCL2","STAT1","MMP2"))

corr_plot <- function(df, title, subtitle, file, w = 9, h = 4.8) {
  lv <- c(intersect(cell_order, df$cell), setdiff(unique(df$cell), cell_order))
  df$cell <- factor(df$cell, levels = lv)
  df$gene <- factor(df$gene, levels = gene_order)
  df$lab  <- paste0(sprintf("%.2f", df$rho), "\n", star(df$FDR))
  p <- ggplot(df, aes(cell, gene, fill = rho)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = lab), size = 3.2, lineheight = 0.9) +
    scale_fill_gradient2(low = "#2E6F9E", mid = "white", high = "#C0392B",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 40, hjust = 1),
          plot.title = element_text(face = "bold"))
  ggsave(file, p, width = w, height = h, dpi = 300)
  print(p)
}

ov <- read.csv("results/hub_immune_corr_overall_FDR.csv")
wi <- read.csv("results/hub_immune_corr_within_FDR.csv")
note <- "BH-adjusted: * < 0.05, ** < 0.01, *** < 0.001"

corr_plot(ov, "Hub genes and immune-cell enrichment: all samples",
          paste("Patients and controls pooled; Spearman correlation;", note),
          "results/hub_immune_corr_v3.png")
corr_plot(wi[wi$group == "UC", ], "Hub genes and immune-cell enrichment: UC patients",
          paste("UC patients only (n = 74); Spearman correlation;", note),
          "results/hub_immune_corr_UC_v3.png")
corr_plot(wi[wi$group == "CD", ], "Hub genes and immune-cell enrichment: CD patients",
          paste("CD patients only (n = 51); Spearman correlation;", note),
          "results/hub_immune_corr_CD_v3.png")