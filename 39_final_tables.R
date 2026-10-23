library(tidyverse)
tab <- read_csv("results/dgidb_gene_drug_summary.csv", show_col_types = FALSE)

# Jadval A: tasdiqlangan biologiklar, hub genlar bilan bevosita bog'liq
A <- tab %>%
  filter(drug_name %in% c("CANAKINUMAB","RILONACEPT","SILTUXIMAB")) %>%
  transmute(Gene = gene_name, Drug = str_to_title(drug_name),
            `Interaction type` = types, `DGIdb sources` = sources,
            `n sources` = n_sources, PMID = "", `Verified (ChEMBL/DrugBank)` = "")

# Jadval B: MMP2 ingibitorlari (docking ligandlari)
lig <- c("MARIMASTAT","BATIMASTAT","ILOMASTAT","PRINOMASTAT","CGS-27023A",
         "ARP100","SB-3CT","REBIMASTAT","TANOMASTAT")
B <- tab %>%
  filter(gene_name == "MMP2", drug_name %in% lig) %>%
  arrange(desc(n_sources)) %>%
  transmute(Drug = str_to_title(drug_name), `Interaction type` = types,
            `DGIdb sources` = sources, `n sources` = n_sources,
            Approved = approved, PMID = "", `Verified (ChEMBL/DrugBank)` = "")

A; B
write.csv(A, "results/Table_A_biologics.csv", row.names = FALSE)
write.csv(B, "results/Table_B_MMP2_inhibitors.csv", row.names = FALSE)