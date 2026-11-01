# Expression-derived cellular proxy estimation and five discovery sensitivity models.
# Analysis-only public implementation; no graphical output is produced.

# Resolve project-relative paths independently of the caller's working directory.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
project_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
setwd(project_root)
if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(BRETIGEA)
  library(yaml)
})

cat("Cell-expression proxy sensitivity models\n")

# Load parameters and seeds
params <- yaml.load_file("config/params.yaml")
seeds <- yaml.load_file("config/seeds.yaml")
set.seed(seeds$global)

dir.create("results/de", recursive = TRUE, showWarnings = FALSE)

# 1. Load inputs
cat("Loading inputs...\n")
vst_mat <- as.matrix(read.csv("data/interim/vst_blind.csv", row.names = 1, check.names = FALSE))
counts_mat <- readRDS("data/interim/counts.rds")
coldata <- read.csv("data/interim/coldata.csv", row.names = 1, stringsAsFactors = FALSE)
gene_annot <- read.csv("data/interim/gene_annotations.csv", stringsAsFactors = FALSE)
rownames(gene_annot) <- as.character(gene_annot$entrez_id)

deg_df <- read.csv("results/de/deg_significant.csv", stringsAsFactors = FALSE)
deg_symbols <- unique(deg_df$symbol[deg_df$symbol != "" & !is.na(deg_df$symbol)])
cat(sprintf("Loaded %d significant DEGs (%d unique symbols) under the retained discovery criterion\n", nrow(deg_df), length(deg_symbols)))

base_res <- read.csv("results/de/deseq2_full.csv", stringsAsFactors = FALSE)
rownames(base_res) <- as.character(base_res$entrez_id)

# 2. Deconvolution using disjoint marker set
cat("Preparing gene-symbol expression matrix for deconvolution...\n")
annot_mapped <- gene_annot[!is.na(gene_annot$symbol) & gene_annot$symbol != "", ]
common_genes <- intersect(rownames(vst_mat), as.character(annot_mapped$entrez_id))

vst_mapped <- vst_mat[common_genes, ]
gene_symbols <- annot_mapped[common_genes, "symbol"]

row_vars <- apply(vst_mapped, 1, var)
split_idx <- split(seq_along(gene_symbols), gene_symbols)
keep_idx <- sapply(split_idx, function(idxs) idxs[which.max(row_vars[idxs])])

vst_sym_mat <- vst_mapped[keep_idx, ]
rownames(vst_sym_mat) <- names(keep_idx)

# Exclude discovery genes and the original fixed nine-symbol exclusion set.
excluded_marker_symbols <- c(deg_symbols, "CCR2", "HTR2C", "FLT3", "ACHE", "MMP8", "NLRP4", "NLRP3", "STAT4", "STAT3")
disjoint_markers <- markers_df_human_brain[!(markers_df_human_brain$markers %in% excluded_marker_symbols), ]
cat(sprintf("BRETIGEA human markers: %d original -> %d disjoint markers (removed %d overlapping test genes)\n",
            nrow(markers_df_human_brain), nrow(disjoint_markers),
            nrow(markers_df_human_brain) - nrow(disjoint_markers)))

# Estimate relative cell-expression proxy scores using findCells with disjoint marker set
cell_props <- findCells(
  inputMat = vst_sym_mat,
  markers = disjoint_markers,
  nMarker = 40,
  method = "SVD",
  scale = TRUE
)

cell_props_df <- as.data.frame(cell_props)
cell_props_df$sample_id <- rownames(cell_props_df)
cell_props_df$diagnosis <- coldata[rownames(cell_props_df), "diagnosis"]
cell_props_df$sex <- coldata[rownames(cell_props_df), "sex"]

write.csv(cell_props_df, "results/de/composition_props.csv", row.names = FALSE)
cat("Saved relative cell-expression proxy scores to results/de/composition_props.csv\n")

# Relative cell-expression proxy comparisons.
cell_types <- c("ast", "end", "mic", "neu", "oli", "opc")

# Descriptive comparisons of relative cell-expression proxy scores
diff_tests <- list()
for (ct in cell_types) {
  v_ad <- cell_props_df[cell_props_df$diagnosis == "AD", ct]
  v_ctrl <- cell_props_df[cell_props_df$diagnosis == "control", ct]
  tt <- t.test(v_ad, v_ctrl)
  diff_tests[[ct]] <- data.frame(
    cell_type = ct,
    mean_ad = mean(v_ad),
    mean_ctrl = mean(v_ctrl),
    diff_ad_minus_ctrl = mean(v_ad) - mean(v_ctrl),
    p_val_ttest = tt$p.value,
    stringsAsFactors = FALSE
  )
}
diff_df <- do.call(rbind, diff_tests)
print(diff_df[, c("cell_type", "diff_ad_minus_ctrl", "p_val_ttest")])

# 3. Model Series M0 to M4 and Collinearity Diagnostics
# M0: ~ sex + diagnosis (base model)
# M1: ~ sex + neu + diagnosis
# M2: ~ sex + mic + diagnosis
# M3: ~ sex + neu + mic + diagnosis (primary model)
# M4: ~ sex + neu + mic + ast + end + opc + diagnosis (oli dropped as reference cell type)

# Helper function to compute VIF
compute_vif <- function(df_preds) {
  vifs <- numeric(ncol(df_preds))
  names(vifs) <- colnames(df_preds)
  for (j in seq_len(ncol(df_preds))) {
    y_var <- df_preds[, j]
    x_vars <- df_preds[, -j, drop = FALSE]
    if (ncol(x_vars) == 0) {
      vifs[j] <- 1.0
    } else {
      fit <- lm(y_var ~ ., data = x_vars)
      r2 <- summary(fit)$r.squared
      vifs[j] <- if (r2 >= 0.9999) 999.0 else 1.0 / (1.0 - r2)
    }
  }
  return(vifs)
}

# Standardize continuous covariates for DESeq2
coldata_models <- coldata
for (ct in cell_types) {
  coldata_models[[ct]] <- scale(cell_props_df[rownames(coldata_models), ct])[, 1]
}
coldata_models$diagnosis <- factor(coldata_models$diagnosis, levels = c("control", "AD"))
coldata_models$sex <- factor(coldata_models$sex, levels = c("F", "M"))

model_specs <- list(
  M0_base = list(formula = ~ sex + diagnosis,
                 preds = c("sex")),
  M1_neu  = list(formula = ~ sex + neu + diagnosis,
                 preds = c("sex", "neu")),
  M2_mic  = list(formula = ~ sex + mic + diagnosis,
                 preds = c("sex", "mic")),
  M3_neu_mic = list(formula = ~ sex + neu + mic + diagnosis,
                    preds = c("sex", "neu", "mic")),
  M4_all_six = list(formula = ~ sex + neu + mic + ast + end + opc + diagnosis, # oli dropped
                    preds = c("sex", "neu", "mic", "ast", "end", "opc"))
)

model_results <- list()
model_meta <- list()

for (m_id in names(model_specs)) {
  cat(sprintf("\n--- Processing Model %s: %s ---\n", m_id, deparse(model_specs[[m_id]]$formula)))
  
  # Diagnostics: condition number and VIF
  preds_cols <- model_specs[[m_id]]$preds
  df_p <- coldata_models[, preds_cols, drop = FALSE]
  # convert factors to numeric dummy
  df_num <- as.data.frame(model.matrix(~ ., data = df_p)[, -1, drop = FALSE])
  vif_vals <- compute_vif(df_num)
  
  full_mm <- model.matrix(model_specs[[m_id]]$formula, data = coldata_models)
  cond_num <- kappa(full_mm, exact = TRUE)
  max_vif <- max(vif_vals)
  
  cat(sprintf("Condition number: %.2f | Max VIF: %.2f (Threshold: 10)\n", cond_num, max_vif))
  print(round(vif_vals, 2))
  
  # Fit or load model
  if (m_id == "M0_base") {
    # Reuse the original sex- and diagnosis-adjusted discovery estimates.
    m_df <- data.frame(
      entrez_id = base_res$entrez_id,
      symbol = base_res$symbol,
      log2fc_shrunk = base_res$log2fc_shrunk,
      log2fc_unshrunk = base_res$log2fc_unshrunk,
      padj = base_res$padj,
      stringsAsFactors = FALSE
    )
  } else {
    cat(sprintf("Fitting DESeq2 for %s...\n", m_id))
    dds_m <- DESeqDataSetFromMatrix(
      countData = counts_mat,
      colData = coldata_models,
      design = model_specs[[m_id]]$formula
    )
    dds_m <- DESeq(dds_m)
    res_m_unshrunk <- results(dds_m, contrast = c("diagnosis", "AD", "control"))
    coef_m <- grep("diagnosis.*AD.*control", resultsNames(dds_m), value = TRUE)
    res_m_shrunk <- lfcShrink(dds_m, coef = coef_m, type = "apeglm", res = res_m_unshrunk)
    
    m_df <- data.frame(
      entrez_id = rownames(counts_mat),
      symbol = gene_annot[rownames(counts_mat), "symbol"],
      log2fc_shrunk = res_m_shrunk$log2FoldChange,
      log2fc_unshrunk = res_m_unshrunk$log2FoldChange,
      padj = res_m_shrunk$padj,
      stringsAsFactors = FALSE
    )
  }
  rownames(m_df) <- as.character(m_df$entrez_id)
  model_results[[m_id]] <- m_df
  
  # Calculate median retention for DEGs in this model
  deg_eids <- as.character(deg_df$entrez_id)
  ret_deg <- m_df[deg_eids, "log2fc_shrunk"] / base_res[deg_eids, "log2fc_shrunk"]
  med_ret <- median(ret_deg, na.rm = TRUE)
  
  # Record model metadata
  model_meta[[m_id]] <- data.frame(
    model_id = m_id,
    formula = deparse(model_specs[[m_id]]$formula),
    condition_number = cond_num,
    max_vif = max_vif,
    vif_summary = paste(names(vif_vals), round(vif_vals, 2), sep = "=", collapse = "; "),
    vif_over_10 = (max_vif > 10),
    n_sig_padj05 = sum(!is.na(m_df$padj) & m_df$padj < 0.05),
    median_retention = med_ret,
    stringsAsFactors = FALSE
  )
}

models_summary_df <- do.call(rbind, model_meta)
write.csv(models_summary_df, "results/de/composition_models.csv", row.names = FALSE)
cat("\nSaved model diagnostics to results/de/composition_models.csv\n")
print(models_summary_df[, c("model_id", "condition_number", "max_vif", "vif_over_10", "n_sig_padj05", "median_retention")])

# Save the full M3 results used by subsequent evidence integration.
write.csv(model_results$M3_neu_mic, "results/de/deseq2_composition_adjusted.csv", row.names = FALSE)
cat("Saved primary adjusted model (M3) full results to results/de/deseq2_composition_adjusted.csv\n")

# 4. Continuous Attenuation & Retention Percentiles across all DEGs
cat("\nCalculating continuous attenuation and retention percentiles across all 156 DEGs...\n")
deg_eids <- as.character(deg_df$entrez_id)

atten_list <- list()
for (m_id in names(model_specs)) {
  m_df <- model_results[[m_id]]
  
  # For all DEGs
  base_lfc <- base_res[deg_eids, "log2fc_shrunk"]
  adj_lfc <- m_df[deg_eids, "log2fc_shrunk"]
  adj_padj <- m_df[deg_eids, "padj"]
  
  retention <- adj_lfc / base_lfc
  sign_agree <- (sign(adj_lfc) == sign(base_lfc))
  
  # Empirical CDF on retention for this model
  retention_ecdf <- ecdf(retention)
  percentiles <- retention_ecdf(retention) * 100
  
  for (i in seq_along(deg_eids)) {
    eid <- deg_eids[i]
    sym <- deg_df$symbol[i]
    atten_list[[paste(m_id, eid, sep = "_")]] <- data.frame(
      model_id = m_id,
      entrez_id = eid,
      symbol = sym,
      base_lfc_shrunk = base_lfc[i],
      adj_lfc_shrunk = adj_lfc[i],
      retention = retention[i],
      retention_percentile = percentiles[i],
      adj_padj = adj_padj[i],
      sign_agreement = sign_agree[i],
      stringsAsFactors = FALSE
    )
  }
}
atten_df <- do.call(rbind, atten_list)
write.csv(atten_df, "results/de/composition_attenuation.csv", row.names = FALSE)
cat("Saved continuous attenuation records to results/de/composition_attenuation.csv\n")
