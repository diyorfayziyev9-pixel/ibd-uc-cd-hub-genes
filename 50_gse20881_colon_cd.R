# 50_gse20881_colon_cd.R
# Hub genes in colonic Crohn's disease: GSE20881 (Agilent GPL1708, two-colour, log ratios).
# Primary: sigmoid colon, inflamed CD vs healthy. Secondary: all colonic segments, patient-level means.

library(GEOquery); library(dplyr); library(pROC)
options(timeout = 600)

g20  <- getGEO("GSE20881", GSEMatrix = TRUE)
e20  <- g20[[1]]
pd20 <- pData(e20)
fd   <- fData(e20)

loc <- pd20$`anatomic location:ch1`; dis <- pd20$`disease:ch1`
inf <- pd20$`inflammation status:ch1`; pat <- pd20$`patient:ch1`
hub6 <- c("IL6","IL1B","CCL2","STAT1","MMP2","PTGS2")

# gene value = mean of probes with the same GENE_SYMBOL (NA ignored)
gene_matrix <- function(genes) sapply(genes, function(g) {
  p <- which(fd$GENE_SYMBOL == g)
  colMeans(exprs(e20)[p, , drop = FALSE], na.rm = TRUE) })

grp <- ifelse(dis == "healthy", "Control",
              ifelse(dis == "crohns disease" & inf == "inflamed", "CD", NA))
is_col <- loc %in% c("ascending colon", "descending colon", "sigmoid colon")

analyse <- function(dd, genes, label) {
  y <- factor(dd$grp, levels = c("Control", "CD"))
  out <- do.call(rbind, lapply(genes, function(g) {
    ok <- !is.na(dd[[g]]); x <- dd[[g]][ok]; yy <- y[ok]
    r  <- roc(yy, x, levels = c("Control", "CD"), direction = "<", quiet = TRUE)
    ci <- ci.auc(r, method = "delong")
    data.frame(set = label, gene = g, n_ctrl = sum(yy == "Control"), n_cd = sum(yy == "CD"),
               median_diff = median(x[yy == "CD"]) - median(x[yy == "Control"]),
               p = wilcox.test(x ~ yy)$p.value,
               AUC = as.numeric(auc(r)), lo = ci[1], hi = ci[3]) }))
  out$FDR <- p.adjust(out$p, method = "BH")
  out
}

ge <- gene_matrix(hub6)

# 1) sigmoid colon
ks    <- loc == "sigmoid colon" & !is.na(grp)
r_sig <- analyse(data.frame(ge[ks, ], grp = grp[ks]), hub6, "sigmoid")

# 2) all colonic segments, patient-level means
kc  <- is_col & !is.na(grp)
dfc <- data.frame(patient = pat[kc], grp = grp[kc], ge[kc, ]) %>%
  group_by(patient, grp) %>%
  summarise(across(all_of(hub6), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
r_col <- analyse(as.data.frame(dfc), hub6, "all_colon_patient_mean")

res20 <- rbind(r_sig, r_col)
write.csv(res20, "results/validation_CD_colon_GSE20881.csv", row.names = FALSE)
print(res20, digits = 3)

# 3) positive controls for platform sensitivity (sigmoid colon)
chk <- c("S100A8","S100A9","LCN2","DUOX2","TIMP1")
chk <- chk[sapply(chk, function(g) any(fd$GENE_SYMBOL == g))]
ge2 <- gene_matrix(chk)
r_chk <- analyse(data.frame(ge2[ks, ], grp = grp[ks]), chk, "sigmoid_positive_controls")
write.csv(r_chk, "results/GSE20881_positive_controls.csv", row.names = FALSE)
print(r_chk, digits = 3)
