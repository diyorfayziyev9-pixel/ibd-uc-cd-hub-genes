common_genes <- read.table("results/common_genes_321.txt", stringsAsFactors = FALSE)$V1

# gen belgisini Entrez ID ga o'tkazish
ids <- bitr(common_genes, fromType = "SYMBOL", toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
nrow(ids)

# fon (universe): platformadagi barcha genlar
bg <- bitr(unique(deg_uc_clean$symbol), fromType = "SYMBOL",
           toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene = ids$ENTREZID, universe = bg$ENTREZID,
                OrgDb = org.Hs.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                qvalueCutoff = 0.05, readable = TRUE)

nrow(as.data.frame(ego))
head(as.data.frame(ego)[, c("Description", "Count", "p.adjust")], 10)

write.csv(as.data.frame(ego), "results/GO_BP_common.csv", row.names = FALSE)