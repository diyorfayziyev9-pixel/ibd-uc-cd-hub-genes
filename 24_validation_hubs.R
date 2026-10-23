library(hgu219.db)
library(AnnotationDbi)

# 1) guruhlarni tekshirish
all(colnames(ex_v) %in% rownames(pdata_v))
grp_raw <- pdata_v$`uc or healthy control:ch1`[match(colnames(ex_v), rownames(pdata_v))]
table(grp_raw)

# 2) probeset -> gen belgisi
sym <- AnnotationDbi::mapIds(hgu219.db, keys = rownames(ex_v),
                             column = "SYMBOL", keytype = "PROBEID",
                             multiVals = "first")
sum(!is.na(sym))

# 3) hub genlar probeset'lari
hub5 <- c("IL6", "IL1B", "CCL2", "STAT1", "MMP2")
probes <- names(sym)[sym %in% hub5]
data.frame(probe = probes, symbol = sym[probes])