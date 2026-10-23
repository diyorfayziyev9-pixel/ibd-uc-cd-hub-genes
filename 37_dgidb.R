library(tidyverse)

dg <- read_tsv("data/dgidb_raw.tsv", show_col_types = FALSE) %>%
  rename(approval = `regulatory approval`, score = `interaction score`) %>%
  mutate(drug = toupper(drug))

# 1) tur bo'yicha avtomatik tasnif
biologic_pat <- "MAB$|CEPT$|KINRA$|TUG$|INTERFERON|RECOMBINANT|PEGOL|VACCINE|BINETRAKIN|THYROGLOBULIN|HEPARAN|INSULIN|LENOGRASTIM|ULINASTATIN"
vague_pat    <- "AGENT$|INHIBITOR$|HERBS|EXTRACT|CREAM|DRUGSATFDA|CHEMBL:|MODULATOR|IMMUNOSUPPRESSANT|CARTILAGE"

dg <- dg %>%
  mutate(class = case_when(
    str_detect(drug, biologic_pat) ~ "Biologic",
    str_detect(drug, vague_pat)    ~ "Unspecified",
    TRUE                           ~ "Small molecule"
  )) %>%
  group_by(gene, drug) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup()

# 2) umumiy xulosa
table(dg$gene, dg$class)

# 3) har gen uchun kichik molekulalar (ball bo'yicha)
small <- dg %>% filter(class == "Small molecule") %>%
  arrange(gene, desc(score))
small %>% group_by(gene) %>% slice_head(n = 8) %>% print(n = 40)

# 4) biologiklar
dg %>% filter(class == "Biologic", approval == "Approved") %>%
  arrange(gene, desc(score)) %>% select(gene, drug, indication, score) %>% print(n = 30)

write.csv(dg, "results/dgidb_classified.csv", row.names = FALSE)
write.csv(small, "results/dgidb_small_molecules.csv", row.names = FALSE)