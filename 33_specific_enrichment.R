library(clusterProfiler); library(org.Hs.eg.db)

uc_only <- setdiff(deg_uc_sig$symbol, deg_cd_sig$symbol)
cd_only <- setdiff(deg_cd_sig$symbol, deg_uc_sig$symbol)
length(uc_only); length(cd_only)

run_enrich <- function(genes, name) {
  ids <- bitr(genes, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
  go <- enrichGO(ids$ENTREZID, universe = bg$ENTREZID, OrgDb = org.Hs.eg.db,
                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05, readable = TRUE)
  write.csv(as.data.frame(go), paste0("results/GO_BP_", name, ".csv"), row.names = FALSE)
  p <- dotplot(go, showCategory = 12) + ggtitle(paste("GO BP:", gsub("_", "-", name)))
  ggsave(paste0("results/GO_BP_", name, ".png"), p, width = 9, height = 7, dpi = 300)
  print(head(as.data.frame(go)[, c("Description", "Count", "p.adjust")], 8))
}

run_enrich(uc_only, "UC_specific")
run_enrich(cd_only, "CD_specific")
