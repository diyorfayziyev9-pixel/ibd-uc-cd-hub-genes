library(limma)

group_uc <- factor(group_uc, levels = c("control", "ulcerative colitis"),
                   labels = c("control", "UC"))

design_uc <- model.matrix(~ 0 + group_uc)
colnames(design_uc) <- levels(group_uc)

fit_uc <- lmFit(ex_uc, design_uc)
contrast_uc <- makeContrasts(UC - control, levels = design_uc)
fit_uc2 <- contrasts.fit(fit_uc, contrast_uc)
fit_uc2 <- eBayes(fit_uc2)

deg_uc <- topTable(fit_uc2, number = Inf, adjust.method = "BH")
head(deg_uc)
dim(deg_uc)
sum(deg_uc$adj.P.Val < 0.05 & abs(deg_uc$logFC) > 1)