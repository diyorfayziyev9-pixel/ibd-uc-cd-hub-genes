# 48_hub_methods_comparison.R
# Hub selection with alternative centrality measures on the largest STRING component.
# Requires: g_big (largest component, score >= 0.7; created in 18_hubs_fix.R)

library(igraph)
hub5 <- c("IL6","IL1B","CCL2","STAT1","MMP2")

d <- data.frame(gene = V(g_big)$name,
                degree      = degree(g_big),
                betweenness = betweenness(g_big),
                closeness   = closeness(g_big),
                eigenvector = eigen_centrality(g_big)$vector,
                pagerank    = page_rank(g_big)$vector)

# MCC (cytoHubba definition): each maximal clique adds (size - 1)! to its nodes
cl  <- max_cliques(g_big)
mcc <- setNames(rep(0, vcount(g_big)), V(g_big)$name)
for (cq in cl) { s <- factorial(length(cq) - 1)
                 for (v in as_ids(cq)) mcc[v] <- mcc[v] + s }
d$MCC <- mcc[d$gene]

meths <- c("degree","betweenness","closeness","eigenvector","pagerank","MCC")
rk <- sapply(meths, function(m) rank(-d[[m]], ties.method = "min"))
rownames(rk) <- d$gene
d$mean_rank <- rowMeans(rk)

print(rk[hub5, ])
print(rowSums(rk[hub5, ] <= 10))
print(head(d[order(d$mean_rank), c("gene", "mean_rank")], 10), row.names = FALSE)

# hubs under inclusive handling of ties at the 10th position
t3 <- lapply(c("degree","betweenness","closeness"),
             function(m) d$gene[rank(-d[[m]], ties.method = "min") <= 10])
print(Reduce(intersect, t3))

out <- d; out[paste0("rank_", meths)] <- rk
write.csv(out, "results/hub_methods_comparison.csv", row.names = FALSE)
