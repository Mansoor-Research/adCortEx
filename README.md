# Cortical gene expression in Alzheimer's disease

Analysis code for **Presynaptic and neuropeptide gene expression shows concordance and cellular heterogeneity in Alzheimer's disease cortex**.

This repository reproduces the statistical tables underlying the study: cortical differential expression, the fixed protein-coding signature, external-cohort agreement, age and cellular-proxy sensitivities, Reactome enrichment, reference-cell expression, published cell-population pathology associations and descriptive evidence integration. It contains no source datasets, participant metadata, fitted results, manuscript-production code or figure-generation scripts.

## Software environment

Use R 4.4.3 with Bioconductor 3.20 and Python 3.12. The retained R package versions and their dependency closure are recorded in `renv.lock`; Python dependencies are in `requirements.txt`. From the repository root:

```sh
python -m venv .venv
# Activate .venv using the command appropriate to your operating system.
python -m pip install -r requirements.txt
Rscript -e "install.packages('renv', repos='https://cloud.r-project.org'); renv::restore(project='.', prompt=FALSE)"
```

Add `Rscript` to PATH, or pass its executable path with `--rscript`. Analysis scripts locate the repository from their own file paths rather than requiring a particular machine or user directory. R scripts activate the project library when `renv/activate.R` exists. The R environment check stops on a mismatch in the R version or a direct analysis package.

## Inputs and source versions

`config/source_manifest.json` records the input URLs, file sizes and SHA-256 checksums. The initial inputs are the GEO-supplied GSE125583 count archive and series matrix, the GSE33000 series matrix, the GPL4372 annotation, and five Reactome/Human Protein Atlas annotation files. Discovery metadata are parsed from the supplied series matrix; there is no read alignment. The published ROSMAP `gpath` tables are retrieved from `mathyslab7/ROSMAP_snRNAseq_PFC` at commit `28388f829be6affa1578e56fa09c4be8cd273817`. All 56 tables are required; source Git blob identities and SHA-256 checksums are verified.

GEO sources: [GSE125583](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125583), [GSE33000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33000), [GPL4372](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL4372). Annotation sources: [Reactome](https://reactome.org/download-data), [Human Protein Atlas](https://www.proteinatlas.org/about/download). Published population results: [Mathys laboratory repository](https://github.com/mathyslab7/ROSMAP_snRNAseq_PFC).

Some annotation URLs expose a changing current release. The script accepts a file only when it matches the recorded analysis snapshot. If the publisher has replaced it, obtain the exact archived snapshot from the corresponding author or an appropriate source archive. Substituting newer annotations would constitute a different analysis. A mismatch stops the workflow; the code does not silently update the source or select another endpoint.

To retrieve the recorded inputs and run the workflow:

```sh
python run_analysis.py --fetch
```

Alternatively, import an existing collection of the exact source snapshots. The input directory must contain the relative paths listed in the manifest, beginning with `data/`:

```sh
python run_analysis.py --input-dir /path/to/input-snapshots
```

To validate an already prepared environment and the nine initial source files without fitting models:

```sh
python run_analysis.py --check
```

Published population tables are checked by the separate `cell_fetch` and `cell_evidence` stages. The preflight check does not claim that those stages or model fits have completed. A downloaded ZIP of this repository can be used without Git; unavailable repository revision information is recorded as null.

## Execution order and outputs

The runner executes the following stages in dependency order. Outputs are created under ignored `data/`, `results/` and `env/` directories. Existing research outputs should be kept outside this checkout; rerunning a stage writes its analysis outputs again.

| Stage | Implementation | Principal output |
|---|---|---|
| Discovery inputs | `R/01_data_acquisition.R` | Parsed GEO metadata and extracted counts |
| Discovery QC | `R/02_qc_covariates.R` | Filtered count matrix, covariates, PCA tables and QC metrics |
| Differential expression | `R/03_differential_expression.R` | Genome-wide DESeq2 results and selected discovery genes |
| Protein-coding signature | `R/04_protein_coding_filter.R` | `results/de/deg_protein_coding.csv` |
| Cellular proxy models | `R/03b_cell_composition.R` | Relative proxy scores, model diagnostics and effect-retention tables |
| External metadata | `py/11_prepare_gse33000_metadata.py` | Included sample metadata and probe mappings |
| External comparison | `R/11_external_validation.R` | Age/sex-adjusted limma results and fixed-gene agreement |
| External sensitivity | `R/11_external_overlap_sensitivity.R` | Overlap-weighted estimates and covariate balance |
| Biological annotation | `py/12_biological_evidence.py` | Reactome tests, reference expression and exploratory partner tables |
| Annotation verification | `py/12_verify_biological_evidence.py` | Source and enrichment checks |
| Published population retrieval | `py/21_cell_evidence.py fetch` | Pinned source tables and provenance |
| Published population comparison | `py/21_cell_evidence.py analyze` | Complete 118-by-56 result grid and evidence table |
| Age sensitivity | `R/12_discovery_age_sensitivity.R` | Three matched complete-case DESeq2 models |
| Evidence integration | `py/23_verify_age_and_integrate.py` | Verified combined evidence table and five descriptive categories |
| Reported-count verification | `py/verify_results.py` | Checks of the principal reported denominators and counts |

Individual stages may be run after their upstream outputs exist, for example:

```sh
python run_analysis.py --steps cell_evidence integration verify_results
```

The integration stage also requires the previously generated age-model outputs and input-hash record. `--steps` is an execution convenience; it does not create or replace missing upstream results.

## Retained analytical definitions

The discovery filter is counts of at least 10 in at least 70 samples. DESeq2 uses `~ sex + diagnosis`, with female and control reference levels. Membership requires BH-adjusted p < 0.05 and absolute **unshrunk** log2 fold change >= 1. apeglm shrinkage describes effect sizes; it does not define membership. DESeq2 `results()` defaults, including the independent-filter optimization alpha of 0.10, are retained. The discovery threshold was selected retrospectively during study reconstruction; this repository does not characterize the overall study as prospectively preregistered.

The original 156 discovery genes and the fixed additional symbols ACHE, CCR2, FLT3, HTR2C, MMP8, NLRP3, NLRP4, STAT3 and STAT4 are excluded from the BRETIGEA marker pool. The resulting scores are expression-derived proxies, not measured cell fractions. Five full-cohort sensitivity models are preserved, including omission of the oligodendrocyte score in the broad model. The complete-case analysis uses the same 218 AD and 70 control samples in all three models, excludes one sample with missing age without imputation, and does not reselect genes.

External models use the supplied processed array scale, unambiguous probe-to-Entrez mappings and median aggregation within gene. The fixed CCR2 identifier correction from historical Entrez ID 1231 to 729230 is documented in `config/identifier_crosswalk.json`; both relevant probes are retained. External support requires discovery-compatible direction and BH-adjusted p < 0.05. Overlap weighting is a sensitivity analysis of the same external cohort.

Reactome tests use the common tested-gene universe restricted to mapped genes, all eligible pathways of size 10-500, and zero-hit p values of 1 before BH correction. The reference atlas provides expression context rather than AD differential expression. Reactome pathway/complex-derived partner pairs do not establish direct binding. Published `gpath` estimates are reused without raw-nucleus processing or fitting new donor-level models. BH and Bonferroni corrections use all 2,214 finite fixed-gene/population tests; unavailable results remain unavailable. The resulting pair-level support does not establish gene-level FDR control, causality or therapeutic direction.

## Expected checks

With matching source snapshots and software, the reported checks are 118 fixed protein-coding genes; 112 externally mapped genes, 108 direction agreements and 99 direction-compatible BH-supported genes; 102 supported genes in the overlap-weighted sensitivity; and 118, 117 and 76 supported genes across the three matched age models. Published population matching retains 6,608 positions, including 2,214 finite tests. The five integrated categories contain 29, 11, 34, 24 and 20 genes. These are verification expectations, not parameters used to tune the fits or alter thresholds.

Plotting, earlier target-ranking or chemical modeling stages, reference-library production and manuscript assembly are intentionally outside this code release. No credentials or server access are required by the code. Upstream datasets remain subject to their own access and reuse terms.

## Validation of this code release

R and Python syntax checks, the R environment check, and integrity checks for the nine initial input snapshots passed. In a separate checkout, metadata extraction, biological annotation, published population matching and evidence integration were rerun using the archived model estimates. Nine complete regenerated tables matched their retained counterparts within numerical tolerance, including the 6,608-position population table and 118-gene combined evidence table. Principal reported-count checks also passed. Statistical models were not refitted during this packaging exercise; a fresh environment installation and all publisher downloads have not been tested end to end.
