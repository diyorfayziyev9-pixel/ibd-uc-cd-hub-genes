library(GEOquery)
options(timeout = 600)

gse_c  <- getGEO("GSE66407", GSEMatrix = TRUE)
length(gse_c)
eset_c <- gse_c[[1]]

dim(exprs(eset_c))
annotation(eset_c)
pdata_c <- pData(eset_c)
colnames(pdata_c)