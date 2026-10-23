library(GEOquery)
options(timeout = 600)

gse_v  <- getGEO("GSE206171", GSEMatrix = TRUE)
length(gse_v)
eset_v <- gse_v[[1]]

dim(exprs(eset_v))
annotation(eset_v)
pdata_v <- pData(eset_v)
colnames(pdata_v)