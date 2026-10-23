library(tidyverse)
hub5 <- c("IL6","IL1B","CCL2","STAT1","MMP2")

# 1) asosiy dataset (GSE75214): logFC va FDR
uc <- deg_uc_clean %>% filter(symbol %in% hub5) %>%
  transmute(Gene = symbol, `logFC UC` = round(logFC, 2), `FDR UC` = signif(adj.P.Val, 2))
cd <- deg_cd_clean %>% filter(symbol %in% hub5) %>%
  transmute(Gene = symbol, `logFC CD` = round(logFC, 2), `FDR CD` = signif(adj.P.Val, 2))

# 2) mustaqil validatsiya
val <- as_tibble(res) %>%
  mutate(d = ifelse(str_detect(dataset, "UC"), "UC", "CD")) %>%
  transmute(Gene = gene, d, txt = sprintf("FDR %.1e; AUC %.2f", FDR, AUC)) %>%
  pivot_wider(names_from = d, values_from = txt, names_prefix = "Validation ")

# 3) immun hujayralar: eng kuchli bog'langan hujayra va M2 bilan korrelyatsiya
imm <- tibble(
  Gene = hub5,
  `Top correlated cell (Spearman r)` = sapply(hub5, function(g) {
    i <- which.max(cor_m[g, ])
    sprintf("%s (%.2f)", colnames(cor_m)[i], cor_m[g, i]) }),
  `r with M2 macrophages` = round(cor_m[hub5, "Macrophages M2"], 2)
)

# 4) dori holati (Jadval A va B dan; PMID tekshiruvidan keyin yakunlanadi)
drug <- tibble(
  Gene = hub5,
  `Drug evidence (DGIdb)` = c(
    "Siltuximab (approved biologic)",
    "Canakinumab, rilonacept (approved biologics)",
    "No robust small-molecule or approved agent",
    "No robust small-molecule or approved agent",
    "Investigational MMP inhibitors (none approved)")
)

summary_tab <- uc %>% left_join(cd, by = "Gene") %>%
  left_join(val, by = "Gene") %>% left_join(imm, by = "Gene") %>%
  left_join(drug, by = "Gene") %>%
  arrange(match(Gene, hub5))

print(summary_tab, width = Inf)
write.csv(summary_tab, "results/Table_hub_genes_summary.csv", row.names = FALSE)