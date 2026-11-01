# Check the R version and direct analysis dependencies against the retained lock.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
project_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
setwd(project_root)
if (file.exists("renv/activate.R")) source("renv/activate.R")
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install the locked R environment with renv::restore() before running the analysis.")
}
lock <- jsonlite::fromJSON("renv.lock", simplifyVector = FALSE)
stopifnot(as.character(getRversion()) == lock$R$Version)
direct <- c("GEOquery", "digest", "yaml", "DESeq2", "apeglm", "BRETIGEA",
            "org.Hs.eg.db", "AnnotationDbi", "limma", "data.table", "jsonlite")
for (name in direct) {
  if (!requireNamespace(name, quietly = TRUE)) stop("Missing R package: ", name)
  observed <- as.character(packageVersion(name))
  expected <- lock$Packages[[name]]$Version
  if (observed != expected) stop(name, ": expected ", expected, "; installed ", observed)
}
cat("R version and direct package versions match the retained analysis environment.\n")
