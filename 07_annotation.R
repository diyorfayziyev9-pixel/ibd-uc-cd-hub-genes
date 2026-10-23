fdata <- fData(eset)

get_symbol <- function(x) {
  if (is.na(x) || x == "---") return(NA_character_)
  first <- strsplit(x, " /// ", fixed = TRUE)[[1]][1]
  parts <- strsplit(first, " // ", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(NA_character_)
  trimws(parts[2])
}

fdata$symbol <- vapply(fdata$gene_assignment, get_symbol, character(1))

sum(!is.na(fdata$symbol))
head(fdata[!is.na(fdata$symbol), c("ID", "symbol")], 10)