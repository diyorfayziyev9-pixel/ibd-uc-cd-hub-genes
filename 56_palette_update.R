# 56_palette_update.R
# Final colour update for Figures 3, 4, 5, 7 and 8 (more saturated, still readable in greyscale):
#   - heatmaps (Fig 3B,C and Fig 7A): blue - white - red, saturated, values clipped at +/-2.5
#   - GO / KEGG dot plots (Fig 4, Fig 8): sequential dark blue (most significant) to light green
#   - PPI network (Fig 5): darker hub / neighbour colours, darker edges
# Run AFTER 55_figure_fixes.R, in the same R session.
# Objects needed from earlier scripts: ego, kk (scripts 13, 15), deg_uc_sig, deg_cd_sig, bg (script 33),
#   deg_uc_clean, deg_cd_clean, ex_uc, ex_cd, group_uc, group_cd (script 55), avg (script 35).
# Old versions of the changed panels are kept as *_oldpal.png in results/intermediate/.

library(ggplot2); library(pheatmap); library(clusterProfiler); library(org.Hs.eg.db)
library(igraph); library(png); library(patchwork)
source("scripts/43_style_functions.R")   # group_cols()

inter <- "results/intermediate/"
keep_old <- function(f) if (file.exists(f)) file.copy(f, sub("\\.png$", "_oldpal.png", f), overwrite = TRUE)

# ---------------------------------------------------------------- palettes
pal_div <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
brk_div <- seq(-2.5, 2.5, length.out = 101)          # clip extreme values -> more contrast
seq_cols <- c("#08306B", "#2171B5", "#41B6C4", "#A1D99B")   # dark = most significant
# own dot plot (does not depend on enrichplot internals): fill = p.adjust, size = Count
my_dot <- function(res, n, title, wrap = 30) {
  df <- as.data.frame(res)
  df <- head(df[order(df$p.adjust), ], n)
  gr <- strsplit(as.character(df$GeneRatio), "/")
  df$ratio <- sapply(gr, function(x) as.numeric(x[1]) / as.numeric(x[2]))
  df$label <- vapply(df$Description, function(x) paste(strwrap(x, wrap), collapse = "\n"), "")
  df$label <- factor(df$label, levels = unique(df$label[order(df$ratio)]))
  ggplot(df, aes(ratio, label, size = Count, fill = p.adjust)) +
    geom_point(shape = 21, colour = "black") +
    scale_fill_gradientn(colours = seq_cols, trans = "log10", name = "p.adjust",
                         guide = guide_colourbar(reverse = TRUE)) +
    scale_size(range = c(3, 9)) +
    theme_bw(base_size = 17) +
    theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black"),
          plot.title = element_text(size = 17)) +
    labs(x = "GeneRatio", y = NULL, title = title)
}

# ---------------------------------------------------------------- Fig 3B,C heatmaps
new_heatmap2 <- function(deg_clean, ex, groups, file, n = 25) {
  keep_old(file)
  top <- head(deg_clean[order(deg_clean$adj.P.Val), ], n)
  m <- ex[as.character(top$ID), ]
  rownames(m) <- top$symbol
  g <- as.character(groups)
  is_ctrl <- grepl("control|normal|healthy|ctrl", g, ignore.case = TRUE)
  o <- order(!is_ctrl)
  m <- m[, o]; g <- g[o]
  ann <- data.frame(Group = g, row.names = colnames(m))
  png(file, width = 1200, height = 1800, res = 300)
  pheatmap(m, scale = "row", annotation_col = ann,
           annotation_colors = list(Group = group_cols(ann$Group)),
           cluster_cols = FALSE, annotation_names_col = FALSE,
           show_colnames = FALSE, fontsize = 11, fontsize_row = 11,
           treeheight_row = 15, border_color = NA,
           color = pal_div, breaks = brk_div)
  dev.off()
}
new_heatmap2(deg_uc_clean, ex_uc, group_uc, paste0(inter, "heatmap_UC_v2.png"))
new_heatmap2(deg_cd_clean, ex_cd, group_cd, paste0(inter, "heatmap_CD_v2.png"))

# ---------------------------------------------------------------- Fig 4 dot plots
keep_old(paste0(inter, "GO_BP_dotplot.png")); keep_old(paste0(inter, "KEGG_dotplot.png"))
p_go   <- my_dot(ego, 15, "GO Biological Process: UC & CD common DEGs")
p_kegg <- my_dot(kk,  15, "KEGG pathways: UC & CD common DEGs", wrap = 45)
ggsave(paste0(inter, "GO_BP_dotplot.png"),  p_go,   width = 11, height = 8, dpi = 300)
ggsave(paste0(inter, "KEGG_dotplot.png"),   p_kegg, width = 11, height = 8, dpi = 300)

# ---------------------------------------------------------------- Fig 8 dot plots
uc_only <- setdiff(deg_uc_sig$symbol, deg_cd_sig$symbol)
cd_only <- setdiff(deg_cd_sig$symbol, deg_uc_sig$symbol)
spec_plot <- function(genes, name) {
  ids <- bitr(genes, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
  go <- enrichGO(ids$ENTREZID, universe = bg$ENTREZID, OrgDb = org.Hs.eg.db,
                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05, readable = TRUE)
  f <- paste0(inter, "GO_BP_", name, ".png")
  keep_old(f)
  p <- my_dot(go, 12, paste("GO BP:", gsub("_", "-", name)))
  ggsave(f, p, width = 11, height = 8, dpi = 300)
}
spec_plot(uc_only, "UC_specific")
spec_plot(cd_only, "CD_specific")

# ---------------------------------------------------------------- Fig 7A xCell heatmap
m <- avg[, c("Control_colon", "UC", "CD", "Control_ileum")]
colnames(m) <- c("Control (colon)", "UC", "CD", "Control (ileum)")
f <- paste0(inter, "xcell_heatmap_v2.png"); keep_old(f)
graphics.off()
png(f, width = 2400, height = 2400, res = 300)
pheatmap(m, scale = "row", color = pal_div, breaks = brk_div, border_color = NA,
         fontsize = 11, angle_col = 45,
         main = "xCell enrichment scores by group (row-scaled)")
dev.off()

# ---------------------------------------------------------------- Fig 5 network
edges <- read.delim("data/string_interactions_short_283_nodes.tsv", stringsAsFactors = FALSE)
edges <- data.frame(a = edges[[1]], b = edges[[2]])
g <- simplify(graph_from_data_frame(edges, directed = FALSE))
comp <- components(g)
g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
hubs   <- c("IL6", "IL1B", "CCL2", "STAT1", "MMP2")
is_hub <- V(g)$name %in% hubs
nb <- unique(unlist(lapply(hubs, function(h) neighbors(g, h)$name)))
is_nb <- V(g)$name %in% nb & !is_hub
col_hub <- "#B2182B"; col_nb <- "#4C78A8"; col_other <- "#C5CFDB"
V(g)$color <- ifelse(is_hub, col_hub, ifelse(is_nb, col_nb, col_other))
V(g)$size  <- ifelse(is_hub, 15, ifelse(is_nb, 5, 4))
V(g)$frame.color <- ifelse(is_hub, "black", "white")
V(g)$label <- ifelse(is_hub, V(g)$name, NA)
V(g)$label.color <- "white"; V(g)$label.font <- 2; V(g)$label.family <- "sans"
ends_mat  <- ends(g, E(g))
touch_hub <- ends_mat[, 1] %in% hubs | ends_mat[, 2] %in% hubs
E(g)$color <- ifelse(touch_hub, "#8B97A6", "#D5DBE2")
E(g)$width <- ifelse(touch_hub, 0.5, 0.4)
set.seed(7)
n <- vcount(g)
ang <- seq(pi/2, pi/2 + 2*pi*(4/5), length.out = 5)
hub_xy <- cbind(6 * cos(ang), 6 * sin(ang))
rownames(hub_xy) <- c("IL6", "CCL2", "STAT1", "MMP2", "IL1B")
minx <- rep(-Inf, n); maxx <- rep(Inf, n); miny <- rep(-Inf, n); maxy <- rep(Inf, n)
for (h in rownames(hub_xy)) {
  i <- which(V(g)$name == h)
  minx[i] <- maxx[i] <- hub_xy[h, 1]; miny[i] <- maxy[i] <- hub_xy[h, 2]
}
lay <- layout_with_fr(g, niter = 1500, minx = minx, maxx = maxx, miny = miny, maxy = maxy)
f <- paste0(inter, "PPI_network_hubs.png"); keep_old(f)
graphics.off()
png(f, width = 2400, height = 2000, res = 300)
par(mar = c(1, 1, 1, 1))
plot(g, layout = lay, vertex.label.cex = 0.65)
legend("bottomleft", legend = c("Hub gene", "Direct neighbour of a hub", "Other protein"),
       pch = 21, pt.bg = c(col_hub, col_nb, col_other), pt.cex = 2, bty = "n", cex = 0.8)
dev.off()

# ---------------------------------------------------------------- reassemble figures
source("scripts/44_make_figures.R")      # rebuilds Figures 2-8 and S1 from the panels
# 44 uses older panels for Figure 7 and S1; rebuild them with the final panels (as in script 52)
make_fig(c(paste0(inter, "xcell_heatmap_v2.png"), paste0(inter, "hub_immune_corr_v3.png")),
         2, "results/figures/Figure7_immune.png", 16, 7)
make_fig(c(paste0(inter, "hub_immune_corr_UC_v3.png"), paste0(inter, "hub_immune_corr_CD_v3.png")),
         2, "results/figures/FigureS1_within_group_corr.png", 16, 6)
cat("Done. New figures are in results/figures/\n")
