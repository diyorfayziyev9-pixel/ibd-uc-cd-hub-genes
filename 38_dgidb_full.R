library(tidyverse)
ix <- read_tsv("data/interactions.tsv", comment = "#", show_col_types = FALSE)

hub5 <- c("IL6","IL1B","CCL2","STAT1","MMP2")

tab <- ix %>%
  filter(gene_name %in% hub5) %>%
  group_by(gene_name, drug_name) %>%
  summarise(
    approved = any(drug_is_approved),
    types    = paste(sort(unique(na.omit(interaction_types))), collapse = "; "),
    sources  = paste(sort(unique(interaction_source_db_name)), collapse = "; "),
    n_sources = n_distinct(interaction_source_db_name),
    .groups = "drop"
  ) %>%
  mutate(types = ifelse(types == "", "not specified", types),
         drug_name = toupper(drug_name))

# 1) tanlangan docking ligandlari (MMP2)
lig <- c("MARIMASTAT","BATIMASTAT","ILOMASTAT","PRINOMASTAT","CGS-27023A",
         "ARP100","SB-3CT","REBIMASTAT","TANOMASTAT")
tab %>% filter(gene_name == "MMP2", drug_name %in% lig) %>% print(n = 20, width = Inf)

# 2) tasdiqlangan biologiklar
bio <- c("CANAKINUMAB","RILONACEPT","TOCILIZUMAB","SILTUXIMAB","SARILUMAB",
         "INFLIXIMAB","ADALIMUMAB","USTEKINUMAB")
tab %>% filter(drug_name %in% bio) %>% print(n = 30, width = Inf)

write.csv(tab, "results/dgidb_gene_drug_summary.csv", row.names = FALSE)