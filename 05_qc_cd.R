ex_cd    <- ex[, cd_idx]
group_cd <- pdata$disease[cd_idx]

pca_cd <- prcomp(t(ex_cd), scale. = TRUE)
pc_cd  <- as.data.frame(pca_cd$x[, 1:2])
pc_cd$group <- group_cd

ggplot(pc_cd, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA: CD (ileum, active) vs control")