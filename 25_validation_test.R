library(tidyverse)

hub5 <- c("IL6", "IL1B", "CCL2", "STAT1", "MMP2")

# 1) har gen uchun eng yuqori ifodali probeset (AFFX siz)
cand <- names(sym)[sym %in% hub5 & !grepl("^AFFX", names(sym))]
best <- sapply(hub5, function(g) {
  p <- cand[sym[cand] == g]
  p[which.max(rowMeans(ex_v[p, , drop = FALSE]))]
})
best

# 2) namuna jadvali
meta <- data.frame(
  gsm   = colnames(ex_v),
  title = pdata_v$title[match(colnames(ex_v), rownames(pdata_v))],
  grp   = grp_raw
)
meta$patient <- sub("\\.[ab]$", "", meta$title)
head(meta, 6)
length(unique(meta$patient))

# 3) bemor bo'yicha o'rtacha
expr <- t(ex_v[best, ]); colnames(expr) <- hub5
df <- cbind(meta, expr) %>%
  group_by(patient, grp) %>%
  summarise(across(all_of(hub5), mean), .groups = "drop")
table(df$grp)

# 4) UC active vs healthy
d <- df %>% filter(grp %in% c("Healthy control", "UC - Active disease"))
sapply(hub5, function(g) wilcox.test(d[[g]] ~ d$grp)$p.value)

dl <- d %>% pivot_longer(all_of(hub5), names_to = "gene", values_to = "expr")
p <- ggplot(dl, aes(grp, expr, fill = grp)) + geom_boxplot() +
  facet_wrap(~gene, scales = "free_y") + theme_bw() +
  theme(axis.text.x = element_blank()) + labs(x = NULL, y = "log2 expression")
p
ggsave("results/validation_boxplot.png", p, width = 9, height = 6, dpi = 300)