library(ggplot2)

volcano <- function(df, title, file) {
  df$status <- ifelse(df$adj.P.Val < 0.05 & df$logFC > 1, "Up",
                      ifelse(df$adj.P.Val < 0.05 & df$logFC < -1, "Down", "NS"))
  p <- ggplot(df, aes(logFC, -log10(adj.P.Val), colour = status)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_colour_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw() + labs(title = title, x = "log2 fold change", y = "-log10 adj. P")
  print(p)
  ggsave(file, p, width = 7, height = 6, dpi = 300)
}

volcano(deg_uc_clean, "UC (colon, active) vs control", "results/volcano_UC.png")
volcano(deg_cd_clean, "CD (ileum, active) vs control", "results/volcano_CD.png")