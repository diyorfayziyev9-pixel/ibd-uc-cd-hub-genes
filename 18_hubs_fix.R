library(igraph)

comp <- components(g)
big  <- which.max(comp$csize)
g_big <- induced_subgraph(g, which(comp$membership == big))
vcount(g_big)
ecount(g_big)

hubs2 <- data.frame(
  gene        = V(g_big)$name,
  degree      = degree(g_big),
  betweenness = betweenness(g_big),
  closeness   = closeness(g_big)
)

top_deg <- head(hubs2[order(-hubs2$degree), ], 10)
top_bet <- head(hubs2[order(-hubs2$betweenness), ], 10)
top_clo <- head(hubs2[order(-hubs2$closeness), ], 10)

top_deg
top_bet
top_clo

hub_genes <- Reduce(intersect, list(top_deg$gene, top_bet$gene, top_clo$gene))
hub_genes

write.csv(hubs2[order(-hubs2$degree), ], "results/hub_scores_main_component.csv", row.names = FALSE)