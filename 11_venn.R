install.packages("ggVennDiagram")
library(ggVennDiagram)

common_genes <- m$symbol[m$direction != "opposite"]
length(common_genes)
write.table(common_genes, "results/common_genes_321.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

venn_list <- list(UC = deg_uc_sig$symbol, CD = deg_cd_sig$symbol)
p <- ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue")
p

ggsave("results/Venn_UC_CD.png", p, width = 6, height = 5, dpi = 300)