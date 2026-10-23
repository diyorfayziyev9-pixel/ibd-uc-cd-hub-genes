p_go <- dotplot(ego, showCategory = 15) +
  ggtitle("GO Biological Process: UC & CD common DEGs")
p_go

ggsave("results/GO_BP_dotplot.png", p_go, width = 9, height = 7, dpi = 300)