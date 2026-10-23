deg_uc$ID <- rownames(deg_uc)
deg_uc <- merge(deg_uc, fdata[, c("ID", "symbol")], by = "ID", all.x = TRUE)
deg_uc <- deg_uc[order(deg_uc$adj.P.Val), ]

# gen belgisi yo'q probelarni olib tashlash
deg_uc_clean <- deg_uc[!is.na(deg_uc$symbol), ]

# bir gen bir nechta probega ega bo'lsa, eng ahamiyatlisini qoldirish
deg_uc_clean <- deg_uc_clean[!duplicated(deg_uc_clean$symbol), ]

deg_uc_sig <- subset(deg_uc_clean, adj.P.Val < 0.05 & abs(logFC) > 1)

nrow(deg_uc_clean)
nrow(deg_uc_sig)
table(deg_uc_sig$logFC > 0)

write.csv(deg_uc_clean, "results/DEG_UC_all.csv", row.names = FALSE)
write.csv(deg_uc_sig,   "results/DEG_UC_significant.csv", row.names = FALSE)

head(deg_uc_sig[, c("symbol", "logFC", "adj.P.Val")], 10)