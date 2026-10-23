group_cd <- factor(group_cd, levels = c("control", "CD"))

design_cd <- model.matrix(~ 0 + group_cd)
colnames(design_cd) <- levels(group_cd)

fit_cd <- lmFit(ex_cd, design_cd)
contrast_cd <- makeContrasts(CD - control, levels = design_cd)
fit_cd2 <- contrasts.fit(fit_cd, contrast_cd)
fit_cd2 <- eBayes(fit_cd2)

deg_cd <- topTable(fit_cd2, number = Inf, adjust.method = "BH")
deg_cd$ID <- rownames(deg_cd)
deg_cd <- merge(deg_cd, fdata[, c("ID", "symbol")], by = "ID", all.x = TRUE)
deg_cd <- deg_cd[order(deg_cd$adj.P.Val), ]

deg_cd_clean <- deg_cd[!is.na(deg_cd$symbol), ]
deg_cd_clean <- deg_cd_clean[!duplicated(deg_cd_clean$symbol), ]
deg_cd_sig <- subset(deg_cd_clean, adj.P.Val < 0.05 & abs(logFC) > 1)

nrow(deg_cd_clean)
nrow(deg_cd_sig)
table(deg_cd_sig$logFC > 0)

write.csv(deg_cd_clean, "results/DEG_CD_all.csv", row.names = FALSE)
write.csv(deg_cd_sig,   "results/DEG_CD_significant.csv", row.names = FALSE)

head(deg_cd_sig[, c("symbol", "logFC", "adj.P.Val")], 10)