library(org.Hs.eg.db); library(AnnotationDbi); library(tidyverse); library(pROC)

ex_c <- exprs(eset_c)
pc_  <- pdata_c
sel  <- which(pc_$`tissue:ch1` == "ileum" &
                ((pc_$`diagnosis:ch1` == "CD" & pc_$`inflammation:ch1` == "yes") |
                   pc_$`diagnosis:ch1` == "Control"))
length(sel)
grp_c <- ifelse(pc_$`diagnosis:ch1`[sel] == "CD", "CD", "Control")
table(grp_c)

# ENSG -> gen belgisi
ens <- sub("_at$", "", rownames(ex_c))
sym_c <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL",
                keytype = "ENSEMBL", multiVals = "first")
hub5 <- c("IL6", "IL1B", "CCL2", "STAT1", "MMP2")
hub_rows <- sapply(hub5, function(g) {
  r <- which(sym_c == g)
  if (length(r) == 0) NA else r[which.max(rowMeans(ex_c[r, , drop = FALSE]))]
})
hub_rows

dc <- as.data.frame(t(ex_c[hub_rows, sel])); colnames(dc) <- hub5
dc$grp <- grp_c
dc$patient <- pc_$`patient:ch1`[sel]

# bemor bo'yicha o'rtacha
dcp <- dc %>% group_by(patient, grp) %>%
  summarise(across(all_of(hub5), mean), .groups = "drop")
table(dcp$grp)

sapply(hub5, function(g) wilcox.test(dcp[[g]] ~ dcp$grp)$p.value)

dl <- dcp %>% pivot_longer(all_of(hub5), names_to = "gene", values_to = "expr")
p <- ggplot(dl, aes(grp, expr, fill = grp)) + geom_boxplot() +
  facet_wrap(~gene, scales = "free_y") + theme_bw() +
  theme(axis.text.x = element_blank()) + labs(x = NULL, y = "log2 expression")
p
ggsave("results/validation_boxplot_CD.png", p, width = 9, height = 6, dpi = 300)

dcp$y <- ifelse(dcp$grp == "CD", 1, 0)
roc_list <- lapply(hub5, function(g) roc(dcp$y, dcp[[g]], quiet = TRUE))
names(roc_list) <- hub5
sapply(roc_list, auc)
p_roc <- ggroc(roc_list, linewidth = 1) + theme_bw() +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", colour = "grey50") +
  labs(title = "ROC: hub genes, CD ileum vs control (GSE66407)")
p_roc
ggsave("results/validation_ROC_CD.png", p_roc, width = 7, height = 6, dpi = 300)