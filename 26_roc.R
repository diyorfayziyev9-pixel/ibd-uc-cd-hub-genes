
library(pROC)

d$y <- ifelse(d$grp == "UC - Active disease", 1, 0)

roc_list <- lapply(hub5, function(g) roc(d$y, d[[g]], quiet = TRUE))
names(roc_list) <- hub5
sapply(roc_list, auc)
sapply(roc_list, function(r) as.numeric(ci.auc(r))[c(1, 3)])

p_roc <- ggroc(roc_list, linewidth = 1) + theme_bw() +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", colour = "grey50") +
  labs(title = "ROC: hub genes, UC active vs healthy (GSE206171)")
p_roc
ggsave("results/validation_ROC.png", p_roc, width = 7, height = 6, dpi = 300)