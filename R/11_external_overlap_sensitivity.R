# Same-cohort sensitivity using age and sex overlap weights.
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
base <- file.path(root, "results", "external_validation")
out <- file.path(base, "analysis_v2_historical_ids")
meta <- fread(file.path(base, "sample_metadata.csv"), colClasses = "character")
meta <- meta[disease_status %in% c("Alzheimer's disease", "non-demented")]
meta[, diagnosis := factor(ifelse(disease_status == "Alzheimer's disease", "AD", "Control"),
                           levels = c("Control", "AD"))]
meta[, sex := factor(sex)]
meta[, age_years := as.numeric(age_years)]
stopifnot(all(!is.na(meta$age_years)), all(!is.na(meta$sex)))

expr_table <- fread(file.path(out, "gene_expression.csv.gz"),
                    colClasses = list(character = "entrez_id"))
stopifnot(all(meta$gsm %in% names(expr_table)))
sample_cols <- meta$gsm
genes <- as.matrix(expr_table[, ..sample_cols])
storage.mode(genes) <- "double"
rownames(genes) <- expr_table$entrez_id

propensity <- glm(I(diagnosis == "AD") ~ age_years + sex,
                  data = meta, family = binomial())
prob_ad <- as.numeric(predict(propensity, type = "response"))
stopifnot(all(is.finite(prob_ad)), all(prob_ad > 0 & prob_ad < 1))
weights <- ifelse(meta$diagnosis == "AD", 1 - prob_ad, prob_ad)
stopifnot(all(weights > 0))
meta[, `:=`(prob_ad = prob_ad, overlap_weight = weights)]
fwrite(meta[, .(gsm, diagnosis, sex, age_years, prob_ad, overlap_weight)],
       file.path(out, "overlap_sample_weights.csv"))

balance <- rbindlist(lapply(levels(meta$diagnosis), function(group) {
  idx <- meta$diagnosis == group
  w <- weights[idx]
  data.table(
    diagnosis = group,
    n = sum(idx),
    effective_n = sum(w)^2 / sum(w^2),
    age_mean_unweighted = mean(meta$age_years[idx]),
    age_mean_overlap = weighted.mean(meta$age_years[idx], w),
    male_fraction_unweighted = mean(meta$sex[idx] == "male"),
    male_fraction_overlap = weighted.mean(meta$sex[idx] == "male", w)
  )
}))
fwrite(balance, file.path(out, "overlap_covariate_balance.csv"))

design <- model.matrix(~ diagnosis + sex + age_years, data = meta)
stopifnot(qr(design)$rank == ncol(design), "diagnosisAD" %in% colnames(design))
weight_matrix <- matrix(rep(weights, each = nrow(genes)), nrow = nrow(genes))
fit <- eBayes(lmFit(genes, design, weights = weight_matrix))
tt <- as.data.table(topTable(fit, coef = "diagnosisAD", number = Inf,
                            sort.by = "none", adjust.method = "BH"), keep.rownames = "entrez_id")
setnames(tt, c("logFC", "adj.P.Val", "P.Value"),
         c("overlap_ad_control_coef", "overlap_padj", "overlap_pvalue"))
fwrite(tt, file.path(out, "overlap_gene_level_limma.csv"))

discovery <- fread(file.path(root, "results/de/deg_protein_coding.csv"),
                   colClasses = list(character = "entrez_id"))
check <- merge(discovery[, .(entrez_id, symbol, discovery_lfc = log2fc_unshrunk)],
               tt[, .(entrez_id, overlap_ad_control_coef, overlap_pvalue, overlap_padj)],
               by = "entrez_id", all.x = TRUE, sort = FALSE)
check[, mapped := is.finite(overlap_ad_control_coef)]
check[, same_direction := ifelse(mapped,
                                 sign(discovery_lfc) == sign(overlap_ad_control_coef), NA)]
check[, gene_replication := !is.na(same_direction) & same_direction &
        overlap_padj < replication_padj]
fwrite(check, file.path(out, "overlap_gene_concordance.csv"))

usable <- check[mapped & overlap_ad_control_coef != 0]
summary <- data.table(
  mapped_degs = nrow(usable),
  same_direction = sum(usable$same_direction),
  same_direction_fraction = mean(usable$same_direction),
  gene_level_replications_bh = sum(check$gene_replication),
  age_gap_unweighted = diff(balance$age_mean_unweighted),
  age_gap_overlap = diff(balance$age_mean_overlap),
  male_fraction_gap_unweighted = diff(balance$male_fraction_unweighted),
  male_fraction_gap_overlap = diff(balance$male_fraction_overlap),
  min_propensity = min(prob_ad),
  max_propensity = max(prob_ad)
)
fwrite(summary, file.path(out, "overlap_summary.csv"))
capture.output(sessionInfo(), file = file.path(out, "overlap_r_session_info.txt"))
print(balance)
print(summary)
