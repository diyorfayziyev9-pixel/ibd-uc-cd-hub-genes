
library(igraph)

edges <- read.delim("data/string_interactions_short_283_nodes.tsv",
                    stringsAsFactors = FALSE)
colnames(edges)[1:2]

g <- graph_from_data_frame(edges[, 1:2], directed = FALSE)
g <- simplify(g)

vcount(g)
ecount(g)

hubs <- data.frame(
  gene        = V(g)$name,
  degree      = degree(g),
  betweenness = betweenness(g),
  closeness   = closeness(g)
)

top_deg <- head(hubs[order(-hubs$degree), ], 10)
top_bet <- head(hubs[order(-hubs$betweenness), ], 10)
top_clo <- head(hubs[order(-hubs$closeness), ], 10)

top_deg
top_bet
top_clo

intersect(intersect(top_deg$gene, top_bet$gene), top_clo$gene)

write.csv(hubs[order(-hubs$degree), ], "results/hub_scores.csv", row.names = FALSE)