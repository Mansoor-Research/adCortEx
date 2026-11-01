# Discovery differential expression with DESeq2 and apeglm.
# Analysis-only public implementation; no graphical output is produced.

# Resolve project-relative paths independently of the caller's working directory.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
project_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
setwd(project_root)
if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(yaml)
  library(org.Hs.eg.db)
})

cat("=== S3: Differential Expression Analysis ===\n")

# Load parameters and seeds
params <- yaml.load_file("config/params.yaml")
seeds <- yaml.load_file("config/seeds.yaml")
set.seed(seeds$global)

# Output directory
dir.create("results/de", recursive = TRUE, showWarnings = FALSE)

# 1. Load interim inputs
counts_mat <- readRDS("data/interim/counts.rds")
coldata <- read.csv("data/interim/coldata.csv", row.names = 1, stringsAsFactors = FALSE)
gene_annot <- read.csv("data/interim/gene_annotations.csv", stringsAsFactors = FALSE)
rownames(gene_annot) <- as.character(gene_annot$entrez_id)

# Align coldata factor levels
coldata$diagnosis <- factor(coldata$diagnosis, levels = c("control", "AD"))
coldata$sex <- factor(coldata$sex, levels = c("F", "M"))

cat(sprintf("Input counts: %d genes across %d samples\n", nrow(counts_mat), ncol(counts_mat)))
cat(sprintf("Design formula: ~ sex + diagnosis (reference level: %s)\n", levels(coldata$diagnosis)[1]))

# 2. Build DESeqDataSet and run DESeq
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = coldata,
  design = ~ sex + diagnosis
)

cat("Running DESeq(dds)...\n")
dds <- DESeq(dds)

# 3. Extract results for AD vs control
cat("Extracting results(contrast = c('diagnosis', 'AD', 'control'))...\n")
res_unshrunk <- results(dds, contrast = c("diagnosis", "AD", "control"))

# 4. Effect size shrinkage via apeglm
cat("Running apeglm effect size shrinkage...\n")
res_names <- resultsNames(dds)
cat("Available coefficients:", paste(res_names, collapse = ", "), "\n")
diag_coef <- grep("diagnosis.*AD.*control", res_names, value = TRUE)
if (length(diag_coef) == 0) {
  stop("Could not find diagnosis_AD_vs_control in resultsNames(dds): ", paste(res_names, collapse = ", "))
}
cat("Using coefficient for apeglm shrinkage:", diag_coef, "\n")

res_shrunk <- lfcShrink(dds, coef = diag_coef, type = "apeglm", res = res_unshrunk)

# 5. Assemble comprehensive results table
gene_ids <- rownames(counts_mat)
full_res_df <- data.frame(
  entrez_id = gene_ids,
  symbol = gene_annot[gene_ids, "symbol"],
  ensembl_id = gene_annot[gene_ids, "ensembl_id"],
  genename = gene_annot[gene_ids, "genename"],
  baseMean = res_shrunk$baseMean,
  log2fc_unshrunk = res_unshrunk$log2FoldChange,
  lfcSE_unshrunk = res_unshrunk$lfcSE,
  stat = res_unshrunk$stat,
  pvalue = res_shrunk$pvalue,
  padj = res_shrunk$padj,
  log2fc_shrunk = res_shrunk$log2FoldChange,
  lfcSE_shrunk = res_shrunk$lfcSE,
  stringsAsFactors = FALSE
)

# Save full results table
write.csv(full_res_df, "results/de/deseq2_full.csv", row.names = FALSE)
cat("Full DESeq2 results saved to results/de/deseq2_full.csv\n")

# 6. Call significant DEGs (retained discovery criterion)
padj_thresh <- params$de$padj
lfc_thresh <- params$de$abs_log2fc
lfc_target <- params$de$abs_log2fc_applies_to

# The published membership criterion uses the unshrunk diagnosis coefficient.
# Defaults of results(), including independent-filter alpha = 0.10, are retained.
stopifnot(identical(lfc_target, "log2fc_unshrunk"))
deg_mask <- !is.na(full_res_df$padj) & full_res_df$padj < padj_thresh &
            abs(full_res_df$log2fc_unshrunk) >= lfc_thresh

deg_df <- full_res_df[deg_mask, ]
# Order output rows by adjusted p value and absolute shrunken effect size.
deg_df <- deg_df[order(deg_df$padj, -abs(deg_df$log2fc_shrunk)), ]
write.csv(deg_df, "results/de/deg_significant.csv", row.names = FALSE)

n_deg <- nrow(deg_df)
n_up <- sum(deg_df$log2fc_shrunk > 0)
n_down <- sum(deg_df$log2fc_shrunk < 0)
cat(sprintf("Significant DEGs (retained discovery criterion): %d (Up: %d, Down: %d) [Expected 156]\n", n_deg, n_up, n_down))
