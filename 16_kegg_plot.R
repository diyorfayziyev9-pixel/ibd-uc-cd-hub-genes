p_kegg <- dotplot(kk, showCategory = 15) +
  ggtitle("KEGG pathways: UC & CD common DEGs")
p_kegg

ggsave("results/KEGG_dotplot.png", p_kegg, width = 9, height = 7, dpi = 300)