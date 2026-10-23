library(pheatmap)

make_heatmap <- function(deg_clean, ex, ids, groups, file, n = 40) {
  top <- head(deg_clean[order(deg_clean$adj.P.Val), ], n)
  m <- ex[as.character(top$ID), ]
  rownames(m) <- top$symbol
  ann <- data.frame(Group = groups, row.names = colnames(m))
  png(file, width = 2200, height = 2400, res = 300)
  pheatmap(m, scale = "row", annotation_col = ann,
           show_colnames = FALSE, fontsize_row = 7,
           color = colorRampPalette(c("steelblue", "white", "firebrick"))(100))
  dev.off()
}

make_heatmap(deg_uc_clean, ex_uc, uc_idx, as.character(group_uc), "results/heatmap_UC.png")
make_heatmap(deg_cd_clean, ex_cd, cd_idx, as.character(group_cd), "results/heatmap_CD.png")