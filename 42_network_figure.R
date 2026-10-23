library(igraph)

edges <- read.delim("data/string_interactions_short_283_nodes.tsv",
                    stringsAsFactors = FALSE)
edges <- data.frame(a = edges[[1]], b = edges[[2]])

g <- simplify(graph_from_data_frame(edges, directed = FALSE))
comp <- components(g)
g <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
cat("Tugunlar:", vcount(g), " Qirralar:", ecount(g), "\n")

hubs   <- c("IL6", "IL1B", "CCL2", "STAT1", "MMP2")
is_hub <- V(g)$name %in% hubs

# hubga to'g'ridan-to'g'ri bog'langan tugunlar (birinchi qavat qo'shnilar)
nb <- unique(unlist(lapply(hubs, function(h) neighbors(g, h)$name)))
is_nb <- V(g)$name %in% nb & !is_hub

V(g)$color <- ifelse(is_hub, "#C0392B",
                     ifelse(is_nb, "#9FB3CC", "#DDE3EA"))
V(g)$size  <- ifelse(is_hub, 15, 4)
V(g)$frame.color <- ifelse(is_hub, "black", "white")
V(g)$label <- ifelse(is_hub, V(g)$name, NA)
V(g)$label.color  <- "white"
V(g)$label.font   <- 2
V(g)$label.family <- "sans"

ends_mat  <- ends(g, E(g))
touch_hub <- ends_mat[, 1] %in% hubs | ends_mat[, 2] %in% hubs
E(g)$color <- ifelse(touch_hub, "#B5BDC6", "#E8ECF0")
E(g)$width <- 0.4

set.seed(7)
# Hublarni besh burchak shaklida qo'lda joylashtiramiz
set.seed(7)
n <- vcount(g)
ang <- seq(pi/2, pi/2 + 2*pi*(4/5), length.out = 5)
hub_xy <- cbind(6 * cos(ang), 6 * sin(ang))
rownames(hub_xy) <- c("IL6", "CCL2", "STAT1", "MMP2", "IL1B")  # tartib: atrofda aylanma

minx <- rep(-Inf, n); maxx <- rep(Inf, n)
miny <- rep(-Inf, n); maxy <- rep(Inf, n)
for (h in rownames(hub_xy)) {
  i <- which(V(g)$name == h)
  minx[i] <- maxx[i] <- hub_xy[h, 1]
  miny[i] <- maxy[i] <- hub_xy[h, 2]
}

lay <- layout_with_fr(g, niter = 1500,
                      minx = minx, maxx = maxx,
                      miny = miny, maxy = maxy)

png("results/PPI_network_hubs.png", width = 2400, height = 2000, res = 300)
par(mar = c(1, 1, 1, 1))
plot(g, layout = lay, vertex.label.cex = 0.65)
legend("bottomleft",
       legend = c("Hub gene", "Direct neighbour of a hub", "Other protein"),
       pch = 21, pt.bg = c("#C0392B", "#9FB3CC", "#DDE3EA"),
       pt.cex = 2, bty = "n", cex = 0.8)
dev.off()
legend("bottomleft",
       legend = c("Hub gene", "Direct neighbour of a hub", "Other protein"),
       pch = 21, pt.bg = c("#C0392B", "#9FB3CC", "#DDE3EA"),
       pt.cex = 2, bty = "n", cex = 0.8)
dev.off()