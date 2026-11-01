# Discovery expression filtering, covariate resolution and sample quality control.
# Analysis-only public implementation; no graphical output is produced.

# Resolve project-relative paths independently of the caller's working directory.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
project_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
setwd(project_root)
if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(org.Hs.eg.db)
  library(yaml)
})

cat("=== S2: QC, Filtering, and Covariate Resolution ===\n")

# Load parameters and seeds
params <- yaml.load_file("config/params.yaml")
seeds <- yaml.load_file("config/seeds.yaml")
set.seed(seeds$global)

# Output directories
dir.create("data/interim", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

# 1. Build counts matrix from raw counts
cat("Loading raw counts from data/raw/raw_counts/...\n")
raw_files <- list.files("data/raw/raw_counts", pattern = "\\.tsv\\.gz$", full.names = TRUE)
raw_files <- sort(raw_files)

# Extract GSM IDs
sample_ids <- sapply(basename(raw_files), function(x) strsplit(x, "_")[[1]][1])
names(sample_ids) <- NULL

# Read the first file to get gene IDs
first_df <- read.table(gzfile(raw_files[1]), header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
gene_ids <- as.character(first_df$ID_REF)
n_genes <- length(gene_ids)
n_samples <- length(sample_ids)
cat("Detected", n_genes, "genes across", n_samples, "samples.\n")

# Assemble count matrix
counts_mat <- matrix(0L, nrow = n_genes, ncol = n_samples,
                     dimnames = list(gene_ids, sample_ids))

for (i in seq_along(raw_files)) {
  df <- read.table(gzfile(raw_files[i]), header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
  counts_mat[, i] <- as.integer(df$count)
}
cat("Counts matrix assembled. Dimensions:", paste(dim(counts_mat), collapse = " x "), "\n")

# 2. Build coldata and resolve covariates
cat("Loading and processing sample metadata...\n")
meta_raw <- read.csv("data/raw/metadata_raw.csv", row.names = 1, stringsAsFactors = FALSE)

# Match rows of coldata to columns of counts_mat
meta_aligned <- meta_raw[sample_ids, ]

# Parse diagnosis
# Raw diagnosis values: "Alzheimer's disease", "control"
raw_diag <- meta_aligned$diagnosis.ch1
diagnosis <- ifelse(grepl("Alzheimer", raw_diag, ignore.case = TRUE), "AD",
              ifelse(grepl("control", raw_diag, ignore.case = TRUE), "control", NA))

# Parse sex
raw_sex <- meta_aligned$Sex.ch1
# Normalize sex to M/F
annot_sex <- ifelse(raw_sex %in% c("M", "Male"), "M",
              ifelse(raw_sex %in% c("F", "Female"), "F", NA))

# Parse age and braak score
age <- as.numeric(meta_aligned$age.ch1)
braak <- meta_aligned$braak.score.ch1

# 3. Resolve Sex using marker genes
cat("Resolving sex using XIST and Y-chromosome markers...\n")
# Entrez IDs: XIST = 7503; RPS4Y1 = 6192, DDX3Y = 8653, UTY = 7404
xist_counts <- counts_mat["7503", ]
y_counts <- colSums(counts_mat[c("6192", "8653", "7404"), , drop = FALSE])

# Sex inference logic: XIST >> Y indicates Female; Y >> XIST indicates Male
# Normalized log2 ratio: log2(xist + 1) - log2(y + 1)
sex_ratio <- log2(xist_counts + 1) - log2(y_counts + 1)
inferred_sex <- ifelse(sex_ratio > 1.0, "F", ifelse(sex_ratio < -1.0, "M", "Uncertain"))

# Concordance check
has_annot <- !is.na(annot_sex)
concordant <- (annot_sex[has_annot] == inferred_sex[has_annot])
concordance_rate <- mean(concordant)
cat(sprintf("Sex inference concordance with metadata: %.2f%% (%d / %d)\n",
            concordance_rate * 100, sum(concordant), sum(has_annot)))

# For missing sex in metadata (GSM3577605), use inferred sex
resolved_sex <- annot_sex
missing_sex_idx <- which(is.na(resolved_sex))
for (idx in missing_sex_idx) {
  cat(sprintf("Sample %s has missing metadata sex. Inferred: %s (XIST=%d, Y_sum=%d)\n",
              sample_ids[idx], inferred_sex[idx], xist_counts[idx], y_counts[idx]))
  resolved_sex[idx] <- inferred_sex[idx]
}

# Any mismatches reported explicitly
mismatches <- which(has_annot & !concordant)
if (length(mismatches) > 0) {
  cat("Reporting annotated vs inferred sex mismatches:\n")
  for (m in mismatches) {
    cat(sprintf("  Sample %s: Annotated=%s, Inferred=%s (XIST=%d, Y_sum=%d)\n",
                sample_ids[m], annot_sex[m], inferred_sex[m], xist_counts[m], y_counts[m]))
  }
}

# Save sex concordance report
sex_report <- data.frame(
  sample_id = sample_ids,
  annotated_sex = annot_sex,
  inferred_sex = inferred_sex,
  resolved_sex = resolved_sex,
  xist_count = xist_counts,
  y_genes_sum = y_counts,
  log2_ratio = sex_ratio,
  concordant = ifelse(is.na(annot_sex), NA, annot_sex == inferred_sex),
  stringsAsFactors = FALSE
)
write.csv(sex_report, "results/qc/sex_concordance.csv", row.names = FALSE)

# Build clean coldata
# Set diagnosis factor with "control" as reference level!
coldata <- data.frame(
  sample_id = sample_ids,
  title = meta_aligned$title,
  diagnosis = factor(diagnosis, levels = c("control", "AD")),
  sex = factor(resolved_sex, levels = c("F", "M")),
  age = age,
  braak = braak,
  row.names = sample_ids,
  stringsAsFactors = FALSE
)
cat("Diagnosis factor levels:", paste(levels(coldata$diagnosis), collapse = ", "), 
    "(Reference level:", levels(coldata$diagnosis)[1], ")\n")
write.csv(coldata, "data/interim/coldata.csv", row.names = TRUE)

# 4. Expression filtering
# Keep genes with count >= min_count in at least smallest_group_n samples
min_count <- params$expression_filter$min_count
smallest_group_n <- min(table(coldata$diagnosis))
cat(sprintf("Expression filter: count >= %d in >= %d samples (smallest group size)\n",
            min_count, smallest_group_n))

expressed_mask <- rowSums(counts_mat >= min_count) >= smallest_group_n
surviving_genes <- gene_ids[expressed_mask]
cat(sprintf("Expressed genes surviving filter: %d / %d (expected 19,108)\n",
            length(surviving_genes), n_genes))

filtered_counts <- counts_mat[expressed_mask, ]

# Map gene annotations (Entrez -> Symbol, Ensembl, Name)
cat("Mapping gene annotations via org.Hs.eg.db...\n")
gene_annot <- data.frame(
  entrez_id = surviving_genes,
  symbol = mapIds(org.Hs.eg.db, keys = surviving_genes, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"),
  ensembl_id = mapIds(org.Hs.eg.db, keys = surviving_genes, column = "ENSEMBL", keytype = "ENTREZID", multiVals = "first"),
  genename = mapIds(org.Hs.eg.db, keys = surviving_genes, column = "GENENAME", keytype = "ENTREZID", multiVals = "first"),
  stringsAsFactors = FALSE
)
write.csv(gene_annot, "data/interim/gene_annotations.csv", row.names = FALSE)
cat(sprintf("Annotation coverage: %.1f%% have Symbol, %.1f%% have Ensembl ID\n",
            mean(!is.na(gene_annot$symbol)) * 100,
            mean(!is.na(gene_annot$ensembl_id)) * 100))

# Save filtered counts matrix as RDS
saveRDS(filtered_counts, "data/interim/counts.rds")
cat("Filtered counts matrix saved to data/interim/counts.rds\n")

# 5. Quality Control: Library sizes, detected genes, outlier scan
cat("Running sample-level QC metrics...\n")
library_sizes <- colSums(filtered_counts)
detected_genes <- colSums(filtered_counts > 0)

lib_mean <- mean(library_sizes)
lib_sd <- sd(library_sizes)
det_mean <- mean(detected_genes)
det_sd <- sd(detected_genes)

# QC criteria:
# 1) library size not extreme outlier (> 3 SD from mean)
# 2) detected genes not extreme low (> 3 SD below mean)
lib_outlier <- abs(library_sizes - lib_mean) > 3 * lib_sd
det_outlier <- (detected_genes < det_mean - 3 * det_sd)
qc_fail_count <- as.integer(lib_outlier) + as.integer(det_outlier)

qc_metrics <- data.frame(
  sample_id = sample_ids,
  diagnosis = coldata$diagnosis,
  sex = coldata$sex,
  library_size = library_sizes,
  detected_genes = detected_genes,
  lib_outlier_3sd = lib_outlier,
  det_outlier_3sd = det_outlier,
  qc_failures = qc_fail_count,
  qc_pass = (qc_fail_count <= 1),
  stringsAsFactors = FALSE
)
write.csv(qc_metrics, "results/qc/qc_metrics.csv", row.names = FALSE)
cat("Samples failing > 1 QC criterion:", sum(qc_fail_count > 1), "\n")

# 6. Blind VST transformation and PCA
cat("Performing blind VST transformation (vst(blind = TRUE))...\n")
dds_qc <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = coldata,
  design = ~ sex + diagnosis
)

vsd <- vst(dds_qc, blind = TRUE)
vst_mat <- assay(vsd)

# Save vst_blind.csv
write.csv(vst_mat, "data/interim/vst_blind.csv", row.names = TRUE)
cat("Blind VST matrix saved to data/interim/vst_blind.csv\n")

# PCA on top 5,000 most variable genes
var_genes <- order(rowVars(vst_mat), decreasing = TRUE)[1:min(5000, nrow(vst_mat))]
pca_res <- prcomp(t(vst_mat[var_genes, ]), center = TRUE, scale. = FALSE)

var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
pca_var_df <- data.frame(
  PC = paste0("PC", 1:10),
  variance_explained = var_explained[1:10],
  cumulative_variance = cumsum(var_explained[1:10]),
  stringsAsFactors = FALSE
)
write.csv(pca_var_df, "results/qc/pca_variance.csv", row.names = FALSE)
cat(sprintf("PC1 variance: %.2f%%, PC2 variance: %.2f%%\n",
            var_explained[1] * 100, var_explained[2] * 100))

pca_df <- data.frame(
  sample_id = sample_ids,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  PC3 = pca_res$x[, 3],
  PC4 = pca_res$x[, 4],
  diagnosis = coldata$diagnosis,
  sex = coldata$sex,
  age = coldata$age,
  stringsAsFactors = FALSE
)
write.csv(pca_df, "results/qc/pca_scores.csv", row.names = FALSE)

# Test covariate associations with PCs
cat("Testing covariate associations with PC1 and PC2...\n")
assoc_results <- list()
for (pc_name in c("PC1", "PC2", "PC3", "PC4")) {
  pc_val <- pca_df[[pc_name]]
  fit_sex <- lm(pc_val ~ coldata$sex)
  fit_diag <- lm(pc_val ~ coldata$diagnosis)
  fit_age <- lm(pc_val ~ coldata$age)
  
  assoc_results[[pc_name]] <- data.frame(
    PC = pc_name,
    var_explained = var_explained[as.integer(sub("PC", "", pc_name))],
    r2_sex = summary(fit_sex)$r.squared,
    pval_sex = summary(fit_sex)$coefficients[2, 4],
    r2_diagnosis = summary(fit_diag)$r.squared,
    pval_diagnosis = summary(fit_diag)$coefficients[2, 4],
    r2_age = summary(fit_age)$r.squared,
    pval_age = summary(fit_age)$coefficients[2, 4],
    stringsAsFactors = FALSE
  )
}
assoc_df <- do.call(rbind, assoc_results)
write.csv(assoc_df, "results/qc/covariate_pc_associations.csv", row.names = FALSE)
print(assoc_df)

# Check if sex dominates variance
sex_dom <- assoc_df["PC1", "r2_sex"] > assoc_df["PC1", "r2_diagnosis"] ||
           assoc_df["PC2", "r2_sex"] > assoc_df["PC2", "r2_diagnosis"]
cat("Sex associates strongly with major PC variance:", sex_dom, 
    "(Confirms original rationale for design ~ sex + diagnosis)\n")

# Pairwise sample correlation matrix for the retained variable genes.
sample_cor <- cor(vst_mat[var_genes, ], method = "spearman")
write.csv(sample_cor, "results/qc/sample_correlations.csv", row.names = TRUE)

# Validate sample quality, covariate resolution and saved PCA summaries.
cat("\n--- Running Stage S2 Gate Checks ---\n")
gate_qc_pass <- all(qc_metrics$qc_failures <= 1)
gate_sex_pass <- (concordance_rate >= 0.95)
gate_pca_pass <- file.exists("results/qc/pca_variance.csv")

cat("Gate 1 (≤1 QC failure per sample):", gate_qc_pass, 
    sprintf("(Max failures: %d)\n", max(qc_metrics$qc_failures)))
cat("Gate 2 (Sex inference agrees > 95%%):", gate_sex_pass, 
    sprintf("(Actual: %.2f%%)\n", concordance_rate * 100))
cat("Gate 3 (PCA written out):", gate_pca_pass, "\n")

if (!gate_qc_pass || !gate_sex_pass || !gate_pca_pass) {
  stop("Stage S2 Gate Failed!")
}

cat("\n=== S2 ACCEPTANCE GATE: PASSED ===\n")
