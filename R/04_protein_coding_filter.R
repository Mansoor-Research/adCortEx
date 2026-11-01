# Protein-coding membership and stable identifier annotation.
# Analysis-only public implementation; no graphical output is produced.

# Resolve project-relative paths independently of the caller's working directory.
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
project_root <- dirname(dirname(normalizePath(script_file, mustWork = TRUE)))
setwd(project_root)
if (file.exists("renv/activate.R")) source("renv/activate.R")

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(yaml)
})

cat("=== S4: Protein-Coding Filter & Symbol Reconciliation ===\n")

# Load parameters and seeds
params <- yaml.load_file("config/params.yaml")
seeds <- yaml.load_file("config/seeds.yaml")
set.seed(seeds$global)

dir.create("results/de", recursive = TRUE, showWarnings = FALSE)

# Load the retained discovery genes selected at padj < 0.05 and |unshrunk LFC| >= 1.0.
deg_df <- read.csv("results/de/deg_significant.csv", stringsAsFactors = FALSE)
cat(sprintf("Loaded %d significant DEGs from results/de/deg_significant.csv\n", nrow(deg_df)))

# 2. Annotate biotype and filter for protein_coding
cat("Annotating gene biotype via org.Hs.eg.db...\n")
entrez_keys <- as.character(deg_df$entrez_id)

biotype_map <- mapIds(org.Hs.eg.db, keys = entrez_keys, column = "GENETYPE", 
                      keytype = "ENTREZID", multiVals = "first")

deg_df$biotype <- biotype_map[entrez_keys]
cat("Biotype breakdown:\n")
print(table(deg_df$biotype, useNA = "ifany"))

# Filter to protein-coding genes
pc_mask <- (!is.na(deg_df$biotype)) & (deg_df$biotype == "protein-coding")
deg_pc <- deg_df[pc_mask, ]
cat(sprintf("Retained %d protein-coding DEGs out of %d total DEGs\n", nrow(deg_pc), nrow(deg_df)))

# 3. Ensure Ensembl IDs and map UniProt accessions
cat("Mapping Ensembl IDs and UniProt accessions...\n")
pc_eids <- as.character(deg_pc$entrez_id)

# Fetch all UniProt accessions for these genes
uniprot_df <- select(org.Hs.eg.db, keys = pc_eids, columns = c("ENTREZID", "UNIPROT", "ENSEMBL"), keytype = "ENTREZID")

# Retain the accession overrides used in the original annotation mapping.
canonical_uniprot <- c(
  "729230" = "P41597", # CCR2
  "3358"   = "P28335", # HTR2C
  "2322"   = "P36888", # FLT3
  "43"     = "P22303", # ACHE
  "4317"   = "P22894", # MMP8
  "4318"   = "P14780", # MMP9
  "147945" = "Q96MN2", # NLRP4
  "114548" = "Q96P20", # NLRP3
  "6775"   = "P42226", # STAT4
  "6774"   = "P40763"  # STAT3
)

# Preserve the original accession preference rule for identifier annotation.
best_uniprot <- function(eid, avail_uniprots) {
  if (eid %in% names(canonical_uniprot)) {
    return(canonical_uniprot[[eid]])
  }
  valid_ups <- na.omit(unique(avail_uniprots))
  if (length(valid_ups) == 0) return(NA_character_)
  # Prefer a six-character accession pattern; this does not establish review status.
  sp_mask <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$", valid_ups)
  if (any(sp_mask)) {
    return(valid_ups[sp_mask][1])
  }
  return(valid_ups[1])
}

uniprot_split <- split(uniprot_df$UNIPROT, uniprot_df$ENTREZID)
deg_pc$uniprot <- sapply(pc_eids, function(eid) {
  best_uniprot(eid, uniprot_split[[eid]])
})

# Ensure Ensembl IDs are populated
missing_ens <- is.na(deg_pc$ensembl_id) | deg_pc$ensembl_id == ""
if (any(missing_ens)) {
  ens_map <- mapIds(org.Hs.eg.db, keys = pc_eids[missing_ens], column = "ENSEMBL",
                    keytype = "ENTREZID", multiVals = "first")
  deg_pc$ensembl_id[missing_ens] <- ens_map
}

# Coverage checks
n_with_ensembl <- sum(!is.na(deg_pc$ensembl_id) & deg_pc$ensembl_id != "")
n_with_uniprot <- sum(!is.na(deg_pc$uniprot) & deg_pc$uniprot != "")
pct_ensembl <- n_with_ensembl / nrow(deg_pc) * 100
pct_uniprot <- n_with_uniprot / nrow(deg_pc) * 100

cat(sprintf("Ensembl ID coverage: %.2f%% (%d / %d)\n", pct_ensembl, n_with_ensembl, nrow(deg_pc)))
cat(sprintf("UniProt accession coverage: %.2f%% (%d / %d)\n", pct_uniprot, n_with_uniprot, nrow(deg_pc)))

# Format output dataframe
deg_protein_coding <- data.frame(
  symbol = deg_pc$symbol,
  entrez_id = deg_pc$entrez_id,
  ensembl_id = deg_pc$ensembl_id,
  uniprot = deg_pc$uniprot,
  biotype = deg_pc$biotype,
  baseMean = deg_pc$baseMean,
  log2fc_shrunk = deg_pc$log2fc_shrunk,
  log2fc_unshrunk = deg_pc$log2fc_unshrunk,
  padj = deg_pc$padj,
  stat = deg_pc$stat,
  stringsAsFactors = FALSE
)

# Order by padj, then absolute shrunken LFC
deg_protein_coding <- deg_protein_coding[order(deg_protein_coding$padj, -abs(deg_protein_coding$log2fc_shrunk)), ]

write.csv(deg_protein_coding, "results/de/deg_protein_coding.csv", row.names = FALSE)
cat(sprintf("Saved %d protein-coding DEGs to results/de/deg_protein_coding.csv\n", nrow(deg_protein_coding)))
