library(GEOquery)
options(timeout = 600)

gse  <- getGEO("GSE75214", GSEMatrix = TRUE)
eset <- gse[[1]]

dim(exprs(eset))
pdata <- pData(eset)
colnames(pdata)