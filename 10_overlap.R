common <- intersect(deg_uc_sig$symbol, deg_cd_sig$symbol)
length(common)

# yo'nalishi bir xil bo'lganlarni ajratish
m <- merge(deg_uc_sig[, c("symbol", "logFC", "adj.P.Val")],
           deg_cd_sig[, c("symbol", "logFC", "adj.P.Val")],
           by = "symbol", suffixes = c("_UC", "_CD"))

m$direction <- ifelse(m$logFC_UC > 0 & m$logFC_CD > 0, "up_both",
                      ifelse(m$logFC_UC < 0 & m$logFC_CD < 0, "down_both", "opposite"))
table(m$direction)

write.csv(m, "results/common_DEGs_UC_CD.csv", row.names = FALSE)
head(m[order(m$adj.P.Val_UC), ], 10)