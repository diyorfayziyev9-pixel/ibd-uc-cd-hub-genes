library(igraph)
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
library(httr)

hub5 <- c("IL6","IL1B","CCL2","STAT1","MMP2")

# umumiy va bir yo'nalishli genlar, berilgan |log2FC| chegarasi uchun
shared_for <- function(lfc) {
  s_uc <- deg_uc_clean[!is.na(deg_uc_clean$symbol) & deg_uc_clean$adj.P.Val < 0.05 &
                         abs(deg_uc_clean$logFC) > lfc, ]
  s_cd <- deg_cd_clean[!is.na(deg_cd_clean$symbol) & deg_cd_clean$adj.P.Val < 0.05 &
                         abs(deg_cd_clean$logFC) > lfc, ]
  mm <- merge(s_uc[, c("symbol","logFC")], s_cd[, c("symbol","logFC")],
              by = "symbol", suffixes = c("_UC","_CD"))
  mm$symbol[sign(mm$logFC_UC) == sign(mm$logFC_CD)]
}

# STRING API
string_edges <- function(genes, score) {
  r <- POST("https://string-db.org/api/tsv/network",
            body = list(identifiers = paste(genes, collapse = "\r"),
                        species = 9606, required_score = score,
                        caller_identity = "ibd_sensitivity"),
            encode = "form")
  stop_for_status(r)
  df <- read.delim(text = content(r, "text", encoding = "UTF-8"),
                   stringsAsFactors = FALSE)
  unique(df[, c("preferredName_A", "preferredName_B")])
}

# hub genlar: eng katta komponent, degree/betweenness/closeness top 10 kesishmasi
hub_calc <- function(ed) {
  g0  <- simplify(graph_from_data_frame(ed, directed = FALSE))
  cp  <- components(g0); big <- which.max(cp$csize)
  gb  <- induced_subgraph(g0, which(cp$membership == big))
  d   <- data.frame(gene = V(gb)$name, degree = degree(gb),
                    betweenness = betweenness(gb), closeness = closeness(gb))
  t3  <- list(head(d[order(-d$degree), ], 10)$gene,
              head(d[order(-d$betweenness), ], 10)$gene,
              head(d[order(-d$closeness), ], 10)$gene)
  list(n_nodes = vcount(g0), n_lcc = vcount(gb), e_lcc = ecount(gb),
       hubs = Reduce(intersect, t3))
}

# tekshiruv: lfc = 1 da 321 gen chiqishi kerak
print(c(n = length(shared_for(1)), same_as_common_genes = setequal(shared_for(1), common_genes)))

grid <- expand.grid(lfc = c(0.58, 1, 1.5), score = c(400, 700, 900))
out  <- vector("list", nrow(grid))
for (k in seq_len(nrow(grid))) {
  gs <- shared_for(grid$lfc[k])
  ed <- string_edges(gs, grid$score[k]); Sys.sleep(1)
  hc <- hub_calc(ed)
  out[[k]] <- data.frame(lfc = grid$lfc[k], score = grid$score[k],
                         n_shared = length(gs), n_nodes = hc$n_nodes,
                         n_lcc = hc$n_lcc, e_lcc = hc$e_lcc,
                         n_of_5 = sum(hub5 %in% hc$hubs),
                         hubs = paste(hc$hubs, collapse = ","))
  cat("done", k, "of", nrow(grid), "\n")
}
sens <- do.call(rbind, out)
write.csv(sens, "results/sensitivity_grid.csv", row.names = FALSE)
print(sens, row.names = FALSE)