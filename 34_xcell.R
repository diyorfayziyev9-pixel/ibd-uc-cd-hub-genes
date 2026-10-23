library(xCell)

# GSE75214: probeset -> gen belgisi, bir gen = bir qator (o'rtacha)
ex_all <- exprs(eset)
sym_all <- fdata$symbol[match(rownames(ex_all), fdata$ID)]
keep <- !is.na(sym_all)
ex_sym <- rowsum(ex_all[keep, ], group = sym_all[keep]) /
  as.vector(table(sym_all[keep])[unique(sym_all[keep])][match(sort(unique(sym_all[keep])), names(table(sym_all[keep])))])
dim(ex_sym)

# faqat tanlangan namunalar (UC active + CD active + nazorat)
idx <- sort(unique(c(uc_idx, cd_idx)))
xc <- xCellAnalysis(ex_sym[, idx])
dim(xc)
head(rownames(xc), 20)