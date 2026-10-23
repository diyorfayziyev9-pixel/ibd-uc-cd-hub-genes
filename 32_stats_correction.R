library(pROC)

# UC validatsiyasi (GSE206171)
p_uc <- sapply(hub5, function(g) wilcox.test(d[[g]] ~ d$grp)$p.value)
y_uc <- ifelse(d$grp == "UC - Active disease", 1, 0)
auc_uc <- t(sapply(hub5, function(g) {
  r <- roc(y_uc, d[[g]], quiet = TRUE)
  c(AUC = as.numeric(auc(r)), CI_low = ci.auc(r)[1], CI_high = ci.auc(r)[3])
}))

# CD validatsiyasi (GSE66407)
p_cd <- sapply(hub5, function(g) wilcox.test(dcp[[g]] ~ dcp$grp)$p.value)
auc_cd <- t(sapply(hub5, function(g) {
  r <- roc(dcp$y, dcp[[g]], quiet = TRUE)
  c(AUC = as.numeric(auc(r)), CI_low = ci.auc(r)[1], CI_high = ci.auc(r)[3])
}))

res <- rbind(
  data.frame(dataset = "GSE206171 (UC)", gene = hub5, p = p_uc,
             FDR = p.adjust(p_uc, "BH"), round(auc_uc, 3)),
  data.frame(dataset = "GSE66407 (CD)", gene = hub5, p = p_cd,
             FDR = p.adjust(p_cd, "BH"), round(auc_cd, 3))
)
res$p   <- signif(res$p, 3)
res$FDR <- signif(res$FDR, 3)
res
write.csv(res, "results/validation_stats_table.csv", row.names = FALSE)