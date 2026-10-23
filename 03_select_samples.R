library(tidyverse)

pdata$disease  <- pdata$`disease:ch1`
pdata$tissue   <- pdata$`tissue:ch1`
pdata$activity <- pdata$`disease activity:ch1`

uc_idx <- which(pdata$disease == "ulcerative colitis" & pdata$activity == "active" |
                  pdata$disease == "control" & pdata$tissue == "colon")

cd_idx <- which(pdata$disease == "CD" & pdata$activity == "active" |
                  pdata$disease == "control" & pdata$tissue == "ileum")

length(uc_idx)
length(cd_idx)