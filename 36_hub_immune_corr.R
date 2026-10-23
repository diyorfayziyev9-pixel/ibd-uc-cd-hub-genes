hub5 <- c("IL6","IL1B","CCL2","STAT1","MMP2")
h <- ex_sym[hub5, idx]
cells <- c("Plasma cells","B-cells","Neutrophils","Monocytes","Macrophages M1",
           "Macrophages M2","cDC","aDC","Tregs","CD8+ Tcm","CD4+ memory T-cells")
cells <- intersect(cells, rownames(xc))

cor_m <- sapply(cells, function(c) sapply(hub5, function(g)
  cor(h[g, ], xc[c, ], method = "spearman")))
round(cor_m, 2)

png("results/hub_immune_corr.png", width = 2000, height = 1200, res = 300)
pheatmap(cor_m, display_numbers = TRUE, number_format = "%.2f",
         color = colorRampPalette(c("steelblue","white","firebrick"))(100),
         breaks = seq(-1, 1, length.out = 101))
dev.off()