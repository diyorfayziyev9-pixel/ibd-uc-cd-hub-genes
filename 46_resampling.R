library(limma)
set.seed(42)
hub5 <- c("IL6","IL1B","CCL2","STAT1","MMP2")
B <- 200

# nazorat birinchi daraja bo'lishini ta'minlash
fix_grp <- function(g) {
  g <- factor(g)
  ref <- levels(g)[grepl("control|normal|healthy", levels(g), ignore.case = TRUE)][1]
  relevel(g, ref = ref)
}
grp_uc <- fix_grp(group_uc); grp_cd <- fix_grp(group_cd)
print(table(grp_uc)); print(table(grp_cd))

ids_uc <- as.character(deg_uc_clean$ID); ids_cd <- as.character(deg_cd_clean$ID)
sym_uc <- setNames(deg_uc_clean$symbol, ids_uc)
sym_cd <- setNames(deg_cd_clean$symbol, ids_cd)

run_deg <- function(ex, grp, keep_ids, sub = NULL) {
  if (!is.null(sub)) { ex <- ex[, sub, drop = FALSE]; grp <- droplevels(grp[sub]) }
  ex  <- ex[rownames(ex) %in% keep_ids, ]
  des <- model.matrix(~ grp)
  fit <- eBayes(lmFit(ex, des))
  topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
}
sub_idx <- function(grp, frac = 0.8)
  unlist(lapply(split(seq_along(grp), grp), function(i) sample(i, round(frac * length(i)))))
call_sign <- function(tt, symmap) {
  s <- ifelse(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1, sign(tt$logFC), 0)
  names(s) <- symmap[rownames(tt)]
  s
}

# tekshiruv: to'liq namunada asl natijaga yaqinmi (UC ~1281, CD ~498 gen)
full_u <- call_sign(run_deg(ex_uc, grp_uc, ids_uc), sym_uc)
full_c <- call_sign(run_deg(ex_cd, grp_cd, ids_cd), sym_cd)
print(c(UC_full = sum(full_u != 0), CD_full = sum(full_c != 0)))

res_list <- vector("list", B)
for (b in seq_len(B)) {
  tu <- run_deg(ex_uc, grp_uc, ids_uc, sub_idx(grp_uc))
  tc <- run_deg(ex_cd, grp_cd, ids_cd, sub_idx(grp_cd))
  s_u <- call_sign(tu, sym_uc); s_c <- call_sign(tc, sym_cd)
  cs   <- intersect(names(s_u), names(s_c))
  conc <- cs[s_u[cs] != 0 & s_u[cs] == s_c[cs]]
  res_list[[b]] <- list(uc = s_u[hub5], cd = s_c[hub5], n_conc = length(conc),
                        jacc = length(intersect(conc, common_genes)) /
                          length(union(conc, common_genes)))
  if (b %% 20 == 0) cat("resample", b, "of", B, "\n")
}

U <- sapply(res_list, function(x) x$uc[hub5])
C <- sapply(res_list, function(x) x$cd[hub5])
summ <- data.frame(gene = hub5,
                   UC_sig_pct   = round(rowMeans(U != 0) * 100, 1),
                   CD_sig_pct   = round(rowMeans(C != 0) * 100, 1),
                   concordant_both_pct = round(rowMeans(U != 0 & U == C) * 100, 1))
write.csv(summ, "results/resampling_hub_stability.csv", row.names = FALSE)
print(summ, row.names = FALSE)
print(summary(sapply(res_list, function(x) x$n_conc)))
print(summary(sapply(res_list, function(x) x$jacc)))