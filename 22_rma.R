library(oligo)
library(pd.hg.u219)

length(cel_files)

raw <- read.celfiles(cel_files)
raw

eset_rma <- rma(raw)
ex_v <- exprs(eset_rma)

dim(ex_v)
summary(as.vector(ex_v))

# namuna nomlarini GSM raqamiga o'tkazish
colnames(ex_v) <- sub("_.*", "", colnames(ex_v))
head(colnames(ex_v))

saveRDS(ex_v, "data/GSE206171/ex_v_rma.rds")