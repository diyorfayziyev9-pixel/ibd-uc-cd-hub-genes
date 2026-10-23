library(tidyverse); library(pheatmap)

# namuna guruhlari (idx tartibida)
grp_all <- pdata$disease[idx]
tis_all <- pdata$tissue[idx]
lab <- ifelse(grp_all == "ulcerative colitis", "UC",
              ifelse(grp_all == "CD", "CD",
                     ifelse(tis_all == "colon", "Control_colon", "Control_ileum")))
table(lab)

# immun hujayra tiplari (faqat immunga tegishlilari)
immune <- c("B-cells","CD4+ memory T-cells","CD4+ naive T-cells","CD8+ T-cells",
            "CD8+ Tcm","CD8+ Tem","cDC","DC","Macrophages","Macrophages M1",
            "Macrophages M2","Monocytes","Neutrophils","NK cells","Plasma cells",
            "pDC","Tregs","Mast cells","Eosinophils","Basophils","iDC","aDC")
immune <- intersect(immune, rownames(xc))

test_one <- function(case, ctrl) {
  a <- xc[immune, lab == case]; b <- xc[immune, lab == ctrl]
  data.frame(cell = immune,
             p = sapply(immune, function(c) wilcox.test(a[c, ], b[c, ])$p.value),
             diff = rowMeans(a) - rowMeans(b)) %>%
    mutate(FDR = p.adjust(p, "BH")) %>% arrange(FDR)
}
res_uc <- test_one("UC", "Control_colon")
res_cd <- test_one("CD", "Control_ileum")
head(res_uc, 8); head(res_cd, 8)
write.csv(res_uc, "results/xcell_UC.csv", row.names = FALSE)
write.csv(res_cd, "results/xcell_CD.csv", row.names = FALSE)

# heatmap: o'rtacha ball guruhlar bo'yicha
avg <- sapply(unique(lab), function(l) rowMeans(xc[immune, lab == l]))
png("results/xcell_heatmap.png", width = 1800, height = 2000, res = 300)
pheatmap(avg, scale = "row", color = colorRampPalette(c("steelblue","white","firebrick"))(100))
dev.off()