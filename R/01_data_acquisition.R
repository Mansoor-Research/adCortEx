# Parse discovery metadata and extract the GEO-supplied count files.
# No read alignment or count normalization is performed at this stage.
# Resolve project-relative paths independently of the caller's working directory.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
project_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
setwd(project_root)
if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages(library(GEOquery))
raw_dir <- file.path("data", "raw")
archive <- file.path(raw_dir, "GSE125583_RAW.tar")
matrix_file <- file.path(raw_dir, "GSE125583_series_matrix.txt.gz")
stopifnot(file.exists(archive), file.exists(matrix_file))
eset <- getGEO(filename = matrix_file, getGPL = FALSE)
if (is.list(eset)) eset <- eset[[1]]
metadata <- pData(eset)
stopifnot(nrow(metadata) == 289L, !anyDuplicated(rownames(metadata)))
diagnosis <- metadata[["diagnosis:ch1"]]
stopifnot(sum(grepl("Alzheimer", diagnosis, ignore.case = TRUE)) == 219L,
          sum(grepl("control", diagnosis, ignore.case = TRUE)) == 70L)
write.csv(metadata, file.path(raw_dir, "metadata_raw.csv"), row.names = TRUE)

# Validate archive member paths before extracting into the project input directory.
members <- untar(archive, list = TRUE)
members <- chartr("\\", "/", members)
members <- chartr("\\", "/", members)
stopifnot(!any(grepl("^/|^[A-Za-z]:|(^|/)\\.\\.(/|$)", members)))
count_dir <- file.path(raw_dir, "raw_counts")
dir.create(count_dir, recursive = TRUE, showWarnings = FALSE)
untar(archive, exdir = count_dir)
stopifnot(length(list.files(count_dir, pattern = "\\.tsv\\.gz$")) == 289L)
cat("Discovery metadata and count files prepared.\n")
