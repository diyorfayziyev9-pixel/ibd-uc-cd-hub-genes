untar("data/GSE206171/GSE206171_RAW.tar", exdir = "data/GSE206171/CEL")

cel_files <- list.files("data/GSE206171/CEL", pattern = "\\.CEL\\.gz$",
                        full.names = TRUE, ignore.case = TRUE)
length(cel_files)
head(basename(cel_files), 3)