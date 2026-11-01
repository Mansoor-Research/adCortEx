# Matched complete-case age and cellular proxy sensitivity models.
# Analysis-only public implementation; no graphical output is produced.

# Resolve project-relative paths independently of the caller's working directory.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
project_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
setwd(project_root)
if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({library(DESeq2); library(jsonlite); library(yaml)})
out <- file.path('results', 'discovery_age_sensitivity_v1')
cfg <- fromJSON(file.path('config', 'discovery_age_sensitivity_v1.json'))
seeds <- yaml::read_yaml(file.path('config', 'seeds.yaml'))
set.seed(seeds$global)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
sink(file.path(out, 'execution.log'), split=TRUE)
cat('Started UTC:', format(Sys.time(), tz='UTC', usetz=TRUE), '\n')
inputs <- c('data/interim/counts.rds', 'data/interim/coldata.csv', 'results/de/composition_props.csv',
            'results/de/deg_protein_coding.csv', 'config/discovery_age_sensitivity_v1.json',
            'R/12_discovery_age_sensitivity.R')
hashes <- vapply(inputs, function(x) digest::digest(file=x, algo='sha256'), character(1))
write_json(as.list(hashes), file.path(out, 'input_hashes.json'), pretty=TRUE, auto_unbox=TRUE)
session <- file.path('env', 'session_info', format(Sys.time(), '%Y%m%dT%H%M%SZ', tz='UTC'))
dir.create(session, recursive=TRUE, showWarnings=FALSE)
capture.output(sessionInfo(), file=file.path(session, 'discovery_age_R_session.txt'))
capture.output(sessionInfo(), file=file.path(out, 'R_session.txt'))
counts <- readRDS('data/interim/counts.rds')
cd <- read.csv('data/interim/coldata.csv', row.names=1, stringsAsFactors=FALSE)
props <- read.csv('results/de/composition_props.csv', stringsAsFactors=FALSE)
rownames(props) <- props$sample_id
stopifnot(identical(colnames(counts), rownames(cd)), !anyDuplicated(rownames(props)), nrow(counts)==cfg$expected_genes)
cd$age <- as.numeric(cd$age)
excluded <- cd[!is.finite(cd$age), , drop=FALSE]
write.csv(excluded, file.path(out, 'excluded_missing_age.csv'), row.names=FALSE)
cd <- cd[is.finite(cd$age), , drop=FALSE]
cd$diagnosis <- factor(cd$diagnosis, levels=c('control', 'AD'))
cd$sex <- factor(cd$sex, levels=c('F', 'M'))
stopifnot(sum(cd$diagnosis=='AD')==cfg$expected_complete_AD, sum(cd$diagnosis=='control')==cfg$expected_complete_control,
          !anyNA(cd$sex), !anyNA(cd$diagnosis), all(rownames(cd) %in% rownames(props)))
cd$neu <- props[rownames(cd), 'neu']
cd$mic <- props[rownames(cd), 'mic']
cd$age_centered <- cd$age - mean(cd$age)
counts <- counts[, rownames(cd)]
write.csv(cd, file.path(out, 'complete_case_metadata.csv'), row.names=FALSE)
fixed <- read.csv('results/de/deg_protein_coding.csv', colClasses=c(entrez_id='character'))
stopifnot(nrow(fixed)==cfg$expected_fixed_genes, !anyDuplicated(fixed$entrez_id))
diagnostics <- list()
all_fixed <- list()
for (model_id in names(cfg$models)) {
  formula <- as.formula(cfg$models[[model_id]])
  design <- model.matrix(formula, cd)
  stopifnot(qr(design)$rank==ncol(design), all(is.finite(design)))
  z <- design[, -1, drop=FALSE]
  vifs <- sapply(seq_len(ncol(z)), function(j) {
    if (ncol(z)==1) return(1)
    rsq <- summary(lm(z[,j] ~ z[,-j,drop=FALSE]))$r.squared
    1/(1-rsq)
  })
  cat('\nMODEL', model_id, cfg$models[[model_id]], '\n')
  dds <- DESeqDataSetFromMatrix(countData=counts, colData=cd, design=formula)
  dds <- DESeq(dds, quiet=FALSE)
  result <- results(dds, contrast=c('diagnosis','AD','control'))
  tab <- data.frame(entrez_id=rownames(result), as.data.frame(result), stringsAsFactors=FALSE)
  tab$model_id <- model_id
  write.csv(tab, file.path(out, paste0(model_id, '_all_genes.csv')), row.names=FALSE)
  selected <- merge(fixed, tab, by='entrez_id', all.x=TRUE, sort=FALSE, suffixes=c('_original','_sensitivity'))
  stopifnot(nrow(selected)==nrow(fixed))
  selected$direction_agrees <- !is.na(selected$log2FoldChange) & sign(selected$log2FoldChange)==sign(selected$log2fc_unshrunk)
  selected$supported <- selected$direction_agrees & !is.na(selected$padj_sensitivity) & selected$padj_sensitivity<cfg$q_threshold
  all_fixed[[model_id]] <- selected
  diagnostics[[model_id]] <- data.frame(model_id=model_id, formula=cfg$models[[model_id]], n_samples=nrow(cd), n_AD=sum(cd$diagnosis=='AD'),
                                       n_control=sum(cd$diagnosis=='control'), condition_number=kappa(design, exact=TRUE), max_vif=max(vifs),
                                       vif_details=paste(paste(colnames(z), signif(vifs,5), sep='='), collapse='; '),
                                       n_genes=nrow(tab), n_finite_q=sum(is.finite(tab$padj)), fixed_supported=sum(selected$supported),
                                       fixed_direction_agrees=sum(selected$direction_agrees))
  write.csv(mcols(dds)[,c('baseMean','dispGeneEst','dispFit','dispersion','betaConv')], file.path(out,paste0(model_id,'_fit_diagnostics.csv')))
}
write.csv(do.call(rbind, all_fixed), file.path(out, 'fixed_gene_age_sensitivity.csv'), row.names=FALSE)
write.csv(do.call(rbind, diagnostics), file.path(out, 'model_summary.csv'), row.names=FALSE)
cat('\nCompleted UTC:', format(Sys.time(), tz='UTC', usetz=TRUE), '\n')
sink()
