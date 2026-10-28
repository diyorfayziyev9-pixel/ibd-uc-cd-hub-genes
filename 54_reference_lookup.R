# 54_reference_lookup.R
# Looks up bibliographic details for candidate references (Crossref) and the
# publications linked to the GEO series (PubMed). Writes results/intermediate/ref_lookup.txt
library(httr); library(jsonlite)

queries <- c(
  "A framework for oligonucleotide microarray preprocessing Carvalho Irizarry",
  "Evolving gene/transcript definitions significantly alter the interpretation of GeneChip data Dai",
  "Exploration, normalization, and summaries of high density oligonucleotide array probe level data Irizarry",
  "cytoHubba: identifying hub objects and sub-networks from complex interactome Chin",
  "Variability in the analysis of a single neuroimaging dataset by many teams Botvinik-Nezer",
  "Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach DeLong",
  "Thousands of samples are needed to generate a robust gene list for predicting outcome in cancer Ein-Dor",
  "Outcome signature genes in breast cancer: is there a unique set? Ein-Dor",
  "Prediction of cancer outcome with microarrays: a multiple random validation strategy Michiels",
  "The MicroArray Quality Control (MAQC) project shows inter- and intraplatform reproducibility of gene expression measurements Shi",
  "Tackling the widespread and critical impact of batch effects in high-throughput data Leek",
  "Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology Sturm",
  "Gene annotation bias impedes biomedical research Haynes Tomczak Khatri",
  "Blockade of interleukin 6 trans signaling suppresses T-cell resistance against apoptosis in chronic intestinal inflammation Atreya",
  "IL-1 beta mediates chronic intestinal inflammation by promoting the accumulation of IL-17A secreting innate lymphoid cells and CD4+ Th17 cells Coccia",
  "Activation of signal transducer and activator of transcription (STAT) 1 in human chronic inflammatory bowel disease Schreiber",
  "The role of neutrophils during intestinal inflammation Fournier Parkos",
  "Tofacitinib as induction and maintenance therapy for ulcerative colitis Sandborn",
  "Monocyte chemoattractant protein 1 gene expression in intestinal epithelial cells and inflammatory bowel disease mucosa Reinecker",
  "A pilot randomized trial of a human anti-interleukin-6 receptor monoclonal antibody in active Crohn's disease Ito",
  "Matrix metalloproteinase levels are elevated in inflammatory bowel disease Baugh",
  "Role of matrix metalloproteinases in intestinal inflammation Medina Radomski",
  "Mucosal imbalance of interleukin-1 and interleukin-1 receptor antagonist in inflammatory bowel disease Casini-Raggi",
  "IL-6 signaling in inflammatory bowel disease: pathophysiological role and clinical relevance Mudter Neurath",
  "Pediatric Crohn disease patients exhibit specific ileal transcriptome and microbiome signature Haberman",
  "Regional variation in gene expression in the healthy colon is dysregulated in ulcerative colitis Noble",
  "Macrophages in intestinal inflammation and resolution: a potential therapeutic target in IBD Na",
  "Macrophages in intestinal homeostasis and inflammation Bain Mowat",
  "Oncostatin M drives intestinal inflammation and predicts response to tumor necrosis factor-neutralizing therapy in patients with inflammatory bowel disease West",
  "Chemokine expression in IBD: mucosal chemokine expression is related to disease activity and distribution Banks",
  "Cyclooxygenase 2 is induced in colonic epithelial cells in inflammatory bowel disease Singer",
  "Host-microbe interactions have shaped the genetic architecture of inflammatory bowel disease Jostins",
  "Association analyses identify 38 susceptibility loci for inflammatory bowel disease and highlight shared genetic risk across populations Liu",
  "A genome-wide association study identifies IL23R as an inflammatory bowel disease gene Duerr",
  "Infliximab for induction and maintenance therapy for ulcerative colitis Rutgeerts",
  "Maintenance infliximab for Crohn's disease: the ACCENT I randomised trial Hanauer",
  "A functional genomics predictive network model identifies regulators of inflammatory bowel disease Peters"
)

fmt_auth <- function(a) {
  if (is.null(a) || length(a) == 0) return("")
  paste(vapply(a, function(x) {
    g <- if (is.null(x$given)) "" else gsub("[^A-Z]", "", gsub("([a-z])[a-z]*", "", x$given))
    paste0(if (is.null(x$family)) "" else x$family, " ", g)
  }, ""), collapse = ", ")
}
look <- function(q) tryCatch({
  r <- GET("https://api.crossref.org/works",
           query = list(query.bibliographic = q, rows = 1,
                        select = "title,author,container-title,issued,volume,issue,page,article-number,DOI"),
           user_agent("ref-lookup-R"))
  it <- fromJSON(content(r, "text", encoding = "UTF-8"), simplifyVector = FALSE)$message$items[[1]]
  yr <- tryCatch(it$issued$`date-parts`[[1]][[1]], error = function(e) NA)
  pg <- if (!is.null(it$page)) it$page else if (!is.null(it$`article-number`)) it$`article-number` else ""
  paste0("TITLE: ", it[["title"]][[1]], "\n  AUTHORS: ", fmt_auth(it$author),
         "\n  JOURNAL: ", paste(it$`container-title`, collapse = ""), " | ", yr, " | vol ",
         if (is.null(it$volume)) "" else it$volume, " | issue ", if (is.null(it$issue)) "" else it$issue,
         " | pages ", pg, " | DOI ", it$DOI)
}, error = function(e) paste("ERROR:", conditionMessage(e)))

out <- character()
for (i in seq_along(queries)) {
  out <- c(out, paste0("[", i, "] QUERY: ", queries[i]), paste0("  ", look(queries[i])), "")
  Sys.sleep(0.4)
}

# publications linked to the GEO series (PubMed ID from GEO, details from PubMed)
library(GEOquery)
for (gse in c("GSE75214", "GSE206171", "GSE66407", "GSE20881")) {
  pid <- tryCatch(Meta(getGEO(gse, GSEMatrix = FALSE))$pubmed_id, error = function(e) NA)
  out <- c(out, paste0("GEO ", gse, " PubMed ID: ", paste(pid, collapse = ",")))
  if (!all(is.na(pid)) && length(pid) > 0) {
    r <- GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
             query = list(db = "pubmed", id = pid[1], retmode = "json"))
    s <- fromJSON(content(r, "text", encoding = "UTF-8"))$result[[as.character(pid[1])]]
    out <- c(out, paste0("  ", s$title, " | ", paste(s$authors$name, collapse = ", "), " | ", s$source,
                         " ", s$pubdate, "; vol ", s$volume, "(", s$issue, ") ", s$pages))
  }
}
dir.create("results/intermediate", showWarnings = FALSE)
writeLines(out, "results/intermediate/ref_lookup.txt")
file.edit("results/intermediate/ref_lookup.txt")