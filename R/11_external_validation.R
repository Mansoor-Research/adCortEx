# External cortical microarray comparison with age- and sex-adjusted limma models.
# Analysis-only public implementation; no graphical output is produced.

# Resolve project-relative paths independently of the caller's working directory.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
project_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
setwd(project_root)
if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})

root <- normalizePath(".", winslash = "/")
params <- yaml::read_yaml(file.path(root, "config", "params.yaml"))
replication_padj <- params$external_validation$replication_bh_padj
base_out <- file.path(root, "results", "external_validation")
out <- file.path(base_out, "analysis_v2_historical_ids")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

samples <- fread(file.path(base_out, "sample_metadata.csv"), colClasses = "character")
samples <- samples[disease_status %in% c("Alzheimer's disease", "non-demented")]
stopifnot(nrow(samples) == 467L, uniqueN(samples$gsm) == nrow(samples))
samples[, diagnosis := factor(ifelse(disease_status == "Alzheimer's disease", "AD", "Control"),
                              levels = c("Control", "AD"))]
samples[, sex := factor(sex)]
samples[, age_years := as.numeric(age_years)]
stopifnot(all(!is.na(samples$age_years)), all(!is.na(samples$sex)))

mapping <- fread(file.path(base_out, "probe_mapping.csv"), colClasses = "character")
mapping <- mapping[mapping_status == "unique" & !is.na(entrez_id) & entrez_id != ""]
stopifnot(uniqueN(mapping$probe_id) == nrow(mapping))
# The reviewed historical identifier correction is fixed independently of expression outcomes.
crosswalk <- as.data.table(jsonlite::fromJSON(file.path(root, "config/identifier_crosswalk.json")))
stopifnot(nrow(crosswalk) == 1L, crosswalk$historical_entrez_id == "1231",
          crosswalk$current_entrez_id == "729230", crosswalk$symbol == "CCR2")
repair <- mapping$entrez_id == crosswalk$historical_entrez_id &
          mapping$symbol == crosswalk$symbol
stopifnot(sum(repair) == 2L)
mapping[repair, entrez_id := crosswalk$current_entrez_id]
fwrite(mapping, file.path(out, "probe_mapping_resolved.csv"))

source <- file.path(root, "data", "external", "GSE33000", "GSE33000_series_matrix.txt.gz")
raw <- suppressWarnings(fread(source, skip = "!series_matrix_table_begin"))
stopifnot(all(samples$gsm %in% names(raw)), "ID_REF" %in% names(raw))
raw[, ID_REF := as.character(ID_REF)]
matched <- match(mapping$probe_id, raw$ID_REF)
mapping <- mapping[!is.na(matched)]
matched <- matched[!is.na(matched)]
stopifnot(nrow(mapping) > 0L)

sample_cols <- samples$gsm
expr <- as.matrix(raw[matched, ..sample_cols])
storage.mode(expr) <- "double"
rownames(expr) <- mapping$probe_id
groups <- split(seq_len(nrow(mapping)), mapping$entrez_id)
genes <- matrix(NA_real_, nrow = length(groups), ncol = ncol(expr),
                dimnames = list(names(groups), samples$gsm))
for (i in seq_along(groups)) {
  idx <- groups[[i]]
  if (length(idx) == 1L) {
    genes[i, ] <- expr[idx, ]
  } else {
    genes[i, ] <- apply(expr[idx, , drop = FALSE], 2L, median, na.rm = TRUE)
  }
}
genes[is.nan(genes)] <- NA_real_

gene_table <- data.table(entrez_id = rownames(genes))
gene_table <- cbind(gene_table, as.data.table(genes))
fwrite(gene_table, file.path(out, "gene_expression.csv.gz"))

design <- model.matrix(~ diagnosis + sex + age_years, data = samples)
stopifnot(qr(design)$rank == ncol(design), "diagnosisAD" %in% colnames(design))
fit <- eBayes(lmFit(genes, design))
tt <- as.data.table(topTable(fit, coef = "diagnosisAD", number = Inf,
                            sort.by = "none", adjust.method = "BH"), keep.rownames = "entrez_id")
setnames(tt, c("logFC", "adj.P.Val", "P.Value"),
         c("external_ad_control_coef", "external_padj", "external_pvalue"))
tt <- tt[is.finite(external_ad_control_coef) & is.finite(external_pvalue)]
fwrite(tt, file.path(out, "gene_level_limma.csv"))

discovery <- fread(file.path(root, "results", "de", "deg_protein_coding.csv"),
                   colClasses = list(character = "entrez_id"))
stopifnot(nrow(discovery) == 118L)
concordance <- merge(discovery[, .(entrez_id, symbol, discovery_lfc = log2fc_unshrunk,
                                   discovery_padj = padj)],
                     tt[, .(entrez_id, external_ad_control_coef, external_pvalue, external_padj)],
                     by = "entrez_id", all.x = TRUE, sort = FALSE)
concordance[, mapped := is.finite(external_ad_control_coef)]
concordance[, same_direction := ifelse(mapped,
                                      sign(discovery_lfc) == sign(external_ad_control_coef), NA)]
concordance[, gene_replication := !is.na(same_direction) & same_direction &
              external_padj < replication_padj]
fwrite(concordance, file.path(out, "gene_concordance.csv"))

usable <- concordance[mapped & external_ad_control_coef != 0]
sign_test <- binom.test(sum(usable$same_direction), nrow(usable), p = 0.5)
rho <- suppressWarnings(cor.test(usable$discovery_lfc,
                                 usable$external_ad_control_coef,
                                 method = "spearman", exact = FALSE))
values <- as.numeric(expr)
qc <- data.table(
  samples_ad = sum(samples$diagnosis == "AD"),
  samples_control = sum(samples$diagnosis == "Control"),
  probes_matched = nrow(mapping),
  genes_collapsed = nrow(genes),
  genes_estimable = nrow(tt),
  missing_probe_values = sum(is.na(expr)),
  total_probe_values = length(expr),
  expression_q01 = unname(quantile(values, 0.01, na.rm = TRUE)),
  expression_median = median(values, na.rm = TRUE),
  expression_q99 = unname(quantile(values, 0.99, na.rm = TRUE)),
  design_rank = qr(design)$rank,
  discovery_genes_mapped = sum(concordance$mapped),
  discovery_genes_total = nrow(concordance),
  same_direction_count = sum(usable$same_direction),
  direction_denominator = nrow(usable),
  direction_fraction = mean(usable$same_direction),
  direction_ci_low = sign_test$conf.int[1],
  direction_ci_high = sign_test$conf.int[2],
  direction_binomial_p = sign_test$p.value,
  spearman_rho = unname(rho$estimate),
  spearman_p = rho$p.value,
  replicated_genes_bh = sum(concordance$gene_replication)
)
fwrite(qc, file.path(out, "summary.csv"))
capture.output(sessionInfo(), file = file.path(out, "r_session_info.txt"))
print(qc)
