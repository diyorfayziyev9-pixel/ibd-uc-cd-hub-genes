kk <- enrichKEGG(gene = ids$ENTREZID, universe = bg$ENTREZID,
                 organism = "hsa", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

nrow(as.data.frame(kk))
head(as.data.frame(kk)[, c("Description", "Count", "p.adjust")], 10)

write.csv(as.data.frame(kk), "results/KEGG_common.csv", row.names = FALSE)