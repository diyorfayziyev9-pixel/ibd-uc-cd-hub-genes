library(tidyverse)

ex_uc    <- ex[, uc_idx]
group_uc <- pdata$disease[uc_idx]

pca <- prcomp(t(ex_uc), scale. = TRUE)
pc  <- as.data.frame(pca$x[, 1:2])
pc$group <- group_uc

ggplot(pc, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(title = "PCA: UC (colon, active) vs control")