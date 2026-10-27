# 52_figures_final.R
# ORDER: run 44_make_figures.R first, then this script (it overwrites Figure 7 and S1 with the final versions).
# Requires objects from earlier scripts: avg (from 35_xcell_compare.R),
#   ge, ks, grp, res20 (script 50; run it first in the same session).
# Working (non-final) panels go to results/intermediate/; final figures to results/figures/.

library(ggplot2); library(patchwork); library(png); library(pheatmap)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/intermediate", showWarnings = FALSE)

## helper: same definition as in 44_make_figures.R
make_fig <- function(files, ncol, out, w, h) {
  plots <- lapply(files, function(f) {
    img <- readPNG(f)
    a <- nrow(img) / ncol(img)
    ggplot() +
      annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = a, interpolate = TRUE) +
      coord_fixed(xlim = c(0, 1), ylim = c(0, a), expand = FALSE) +
      theme_void() + theme(plot.margin = margin(28, 4, 4, 28))
  })
  fig <- wrap_plots(plots, ncol = ncol) + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 22))
  ggsave(out, fig, width = w, height = h, dpi = 300)
}

## Figure 7A: xCell enrichment by group (row-scaled heatmap)
stopifnot(all(c("Control_colon", "UC", "CD", "Control_ileum") %in% colnames(avg)))
m <- avg[, c("Control_colon", "UC", "CD", "Control_ileum")]
colnames(m) <- c("Control (colon)", "UC", "CD", "Control (ileum)")
pal <- colorRampPalette(c("#2E6F9E", "white", "#C0392B"))(100)
graphics.off()
png("results/intermediate/xcell_heatmap_v2.png", width = 2400, height = 2400, res = 300)
pheatmap(m, scale = "row", color = pal, border_color = NA, fontsize = 11, angle_col = 45,
         main = "xCell enrichment scores by group (row-scaled)")
dev.off()

## Figure 7 and S1 (51_fig_corr_stars.R writes the v3 files into results/)
v3 <- c("hub_immune_corr_v3.png", "hub_immune_corr_UC_v3.png", "hub_immune_corr_CD_v3.png")
for (f in v3) if (file.exists(file.path("results", f)))
  file.rename(file.path("results", f), file.path("results/intermediate", f))

make_fig(c("results/intermediate/xcell_heatmap_v2.png",
           "results/intermediate/hub_immune_corr_v3.png"),
         ncol = 2, out = "results/figures/Figure7_immune.png", w = 16, h = 7)
make_fig(c("results/intermediate/hub_immune_corr_UC_v3.png",
           "results/intermediate/hub_immune_corr_CD_v3.png"),
         ncol = 2, out = "results/figures/FigureS1_within_group_corr.png", w = 16, h = 6)

## Figure S2: sensitivity of hub selection
sg <- read.csv("results/sensitivity_grid.csv", stringsAsFactors = FALSE)
genes <- c("IL6","IL1B","CCL2","STAT1","MMP2","PTGS2","CXCL8")
long <- do.call(rbind, lapply(seq_len(nrow(sg)), function(i) {
  hubs <- strsplit(sg$hubs[i], ",")[[1]]
  data.frame(lfc = sg$lfc[i], score = sg$score[i], gene = genes, hub = genes %in% hubs) }))
freq <- tapply(long$hub, long$gene, sum)[genes]
long$gene <- factor(long$gene, levels = genes, labels = paste0(genes, "\n", freq, "/9"))
long$lfc_lab <- factor(paste0("|log2FC| > ", long$lfc),
                       levels = paste0("|log2FC| > ", c(1.5, 1, 0.58)))
long$score_lab <- factor(paste0("STRING score >= ", long$score / 1000),
                         levels = paste0("STRING score >= ", c(0.4, 0.7, 0.9)))
p <- ggplot(long, aes(gene, lfc_lab, fill = hub)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  facet_wrap(~ score_lab, nrow = 1) +
  scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#E5E8EC"),
                    labels = c("TRUE" = "Hub gene", "FALSE" = "Not a hub"), name = NULL) +
  labs(title = "Sensitivity of hub-gene selection to analysis settings",
       subtitle = "Primary setting: |log2FC| > 1 and STRING score >= 0.7; numbers show in how many of the nine settings a gene was a hub",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"), legend.position = "bottom")
ggsave("results/figures/FigureS2_sensitivity.png", p, width = 13, height = 4.6, dpi = 300)

## Figure S3: hub genes in colonic CD (GSE20881, sigmoid colon)
hub6 <- c("IL6","IL1B","CCL2","STAT1","MMP2","PTGS2")
dd <- data.frame(ge[ks, hub6], grp = grp[ks])
dl <- do.call(rbind, lapply(hub6, function(g)
  data.frame(gene = g, grp = dd$grp, value = dd[[g]])))
dl$grp  <- factor(dl$grp, levels = c("Control", "CD"), labels = c("Healthy", "CD"))
dl$gene <- factor(dl$gene, levels = hub6)
fc  <- intersect(c("FDR6", "FDR"), names(res20))[1]
rs  <- res20[res20$set == "sigmoid", ]
lab <- data.frame(gene = factor(rs$gene, levels = hub6),
                  txt  = paste0("FDR = ", signif(rs[[fc]], 2), "\nAUC = ", sprintf("%.2f", rs$AUC)))
p <- ggplot(dl, aes(grp, value, fill = grp)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, size = 1, alpha = 0.6) +
  geom_text(data = lab, aes(x = 1.5, y = Inf, label = txt),
            inherit.aes = FALSE, vjust = 1.15, size = 3.3) +
  facet_wrap(~ gene, scales = "free_y", nrow = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  scale_fill_manual(values = c("Healthy" = "#5B7FA6", "CD" = "#8E3B6E"), guide = "none") +
  labs(title = "Hub genes in colonic Crohn's disease (GSE20881, sigmoid colon)",
       subtitle = "Inflamed CD (n = 12) vs healthy (n = 27); Wilcoxon test, BH-adjusted across six genes",
       x = NULL, y = "Log ratio (sample vs reference, Agilent)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
ggsave("results/figures/FigureS3_colonic_CD.png", p, width = 10, height = 6.5, dpi = 300)