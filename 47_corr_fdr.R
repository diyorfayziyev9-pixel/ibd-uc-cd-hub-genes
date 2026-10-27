# 47_corr_fdr.R
# Hub gene - immune cell correlations (Spearman) with BH-adjusted P values.
# Requires objects from the earlier scripts: ex_sym, idx, xc (xCell scores), lab (UC/CD/Control labels)

hub5 <- c("IL6","IL1B","CCL2","STAT1","MMP2")
h <- ex_sym[hub5, idx]
cells <- c("Plasma cells","B-cells","Neutrophils","Monocytes","Macrophages M1",
           "Macrophages M2","cDC","aDC","Tregs","CD8+ Tcm","CD4+ memory T-cells")
cells <- intersect(cells, rownames(xc))
stopifnot(identical(colnames(h), colnames(xc)))

# 1) all samples (patients + controls pooled)
res <- expand.grid(gene = hub5, cell = cells, stringsAsFactors = FALSE)
res$rho <- NA_real_; res$p <- NA_real_
for (i in seq_len(nrow(res))) {
  ct <- suppressWarnings(cor.test(as.numeric(h[res$gene[i], ]),
                                  as.numeric(xc[res$cell[i], ]),
                                  method = "spearman", exact = FALSE))
  res$rho[i] <- unname(ct$estimate); res$p[i] <- ct$p.value
}
res$FDR <- p.adjust(res$p, method = "BH")
write.csv(res, "results/hub_immune_corr_overall_FDR.csv", row.names = FALSE)

# 2) within UC and within CD patients (controls excluded)
within_test <- function(lbl) {
  s <- which(lab == lbl)
  out <- expand.grid(gene = hub5, cell = cells, stringsAsFactors = FALSE)
  out$group <- lbl; out$n <- length(s)
  out$rho <- NA_real_; out$p <- NA_real_
  for (i in seq_len(nrow(out))) {
    ct <- suppressWarnings(cor.test(as.numeric(h[out$gene[i], s]),
                                    as.numeric(xc[out$cell[i], s]),
                                    method = "spearman", exact = FALSE))
    out$rho[i] <- unname(ct$estimate); out$p[i] <- ct$p.value
  }
  out$FDR <- p.adjust(out$p, method = "BH")
  out
}
wt <- rbind(within_test("UC"), within_test("CD"))
write.csv(wt, "results/hub_immune_corr_within_FDR.csv", row.names = FALSE)
