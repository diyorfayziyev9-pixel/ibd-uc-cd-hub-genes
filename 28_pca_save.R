pca_plot <- function(pcdf, title, file) {
  p <- ggplot(pcdf, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) + theme_bw() + labs(title = title)
  print(p)
  ggsave(file, p, width = 7, height = 6, dpi = 300)
}

pca_plot(pc,    "PCA: UC (colon, active) vs control", "results/PCA_UC.png")
pca_plot(pc_cd, "PCA: CD (ileum, active) vs control", "results/PCA_CD.png")