# 49_validation_PTGS2.R
# Supplementary validation of PTGS2 (tied for 10th place by degree) next to the five hub genes.
# Requires objects from earlier scripts: ex_v, sym, meta (UC, GSE206171);
#                                       ex_c, sym_c, sel, grp_c, pc_ (CD, GSE66407)

library(dplyr); library(pROC)
hub5 <- c("IL6","IL1B","CCL2","STAT1","MMP2")
hub6 <- c(hub5, "PTGS2")

## ---- UC: GSE206171 ----
cand6 <- names(sym)[sym %in% hub6 & !grepl("^AFFX", names(sym))]
best6 <- sapply(hub6, function(g) {
  p <- cand6[sym[cand6] == g]
  p[which.max(rowMeans(ex_v[p, , drop = FALSE]))] })
expr6 <- t(ex_v[best6, ]); colnames(expr6) <- hub6
df6 <- cbind(meta, expr6) %>%
  group_by(patient, grp) %>%
  summarise(across(all_of(hub6), mean), .groups = "drop")
d6 <- df6 %>% filter(grp %in% c("Healthy control", "UC - Active disease"))

lv <- c("Healthy control", "UC - Active disease")
y  <- factor(d6$grp, levels = lv)
res6 <- do.call(rbind, lapply(hub6, function(g) {
  r  <- roc(y, d6[[g]], levels = lv, direction = "<", quiet = TRUE)
  ci <- ci.auc(r, method = "delong")
  data.frame(gene = g, p = wilcox.test(d6[[g]] ~ y)$p.value,
             AUC = as.numeric(auc(r)), lo = ci[1], hi = ci[3]) }))
res6$FDR6 <- p.adjust(res6$p, method = "BH")
write.csv(res6, "results/validation_UC_with_PTGS2.csv", row.names = FALSE)

## ---- CD ileum: GSE66407 ----
hub_rows6 <- sapply(hub6, function(g) {
  r <- which(sym_c == g)
  if (length(r) == 0) NA else r[which.max(rowMeans(ex_c[r, , drop = FALSE]))] })
dc6 <- as.data.frame(t(ex_c[hub_rows6, sel])); colnames(dc6) <- hub6
dc6$grp <- grp_c; dc6$patient <- pc_$`patient:ch1`[sel]
dcp6 <- dc6 %>% group_by(patient, grp) %>%
  summarise(across(all_of(hub6), mean), .groups = "drop")

lv_c <- c("Control", "CD")
yc   <- factor(dcp6$grp, levels = lv_c)
res6c <- do.call(rbind, lapply(hub6, function(g) {
  r <- roc(yc, dcp6[[g]], levels = lv_c, direction = "<", quiet = TRUE)
  set.seed(42)
  ci <- ci.auc(r, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
  data.frame(gene = g, p = wilcox.test(dcp6[[g]] ~ yc)$p.value,
             AUC = as.numeric(auc(r)), lo = ci[1], hi = ci[3]) }))
res6c$FDR6 <- p.adjust(res6c$p, method = "BH")
write.csv(res6c, "results/validation_CD_with_PTGS2.csv", row.names = FALSE)
