cells <- c("Plasma cells","B-cells","Neutrophils","Monocytes","Macrophages M1",
           "Macrophages M2","cDC","aDC","Tregs","CD8+ Tcm","CD4+ memory T-cells")
cells <- intersect(cells, rownames(xc))
hub5  <- c("IL6","IL1B","CCL2","STAT1","MMP2")
h <- ex_sym[hub5, idx]

within <- function(lbl) {
  s <- which(lab == lbl)
  sapply(cells, function(c) sapply(hub5, function(g)
    cor(h[g, s], xc[c, s], method = "spearman")))
}
cm <- list(UC = within("UC"), CD = within("CD"))
lapply(cm, round, 2)

for (nm in names(cm)) {
  png(paste0("results/hub_immune_corr_", nm, ".png"), width = 2000, height = 1200, res = 300)
  pheatmap::pheatmap(cm[[nm]], display_numbers = TRUE, number_format = "%.2f",
                     color = colorRampPalette(c("steelblue","white","firebrick"))(100),
                     breaks = seq(-1, 1, length.out = 101), main = paste(nm, "only"))
  dev.off()
}