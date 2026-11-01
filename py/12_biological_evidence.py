"""Target-agnostic bulk pathway, healthy-cell and partner evidence tables."""

from __future__ import annotations

import csv
import json
import math
import re
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/external/biological_evidence"
OUT = ROOT / "results/biological_evidence"
OUT.mkdir(exist_ok=True, parents=True)


def clean_id(value: object) -> str:
    return str(value).strip().removesuffix(".0")


def load_genes() -> tuple[pd.DataFrame, set[str], set[str]]:
    discovery = pd.read_csv(ROOT / "results/de/deseq2_full.csv", dtype={"entrez_id": str})
    external = pd.read_csv(
        ROOT / "results/external_validation/analysis_v2_historical_ids/gene_level_limma.csv",
        dtype={"entrez_id": str},
    )
    deg = pd.read_csv(ROOT / "results/de/deg_protein_coding.csv", dtype={"entrez_id": str})
    for df in (discovery, external, deg):
        df["entrez_id"] = df.entrez_id.map(clean_id)
    common = discovery.merge(external, on="entrez_id", how="inner", validate="one_to_one")
    common["replicated_direction"] = (
        common.log2fc_unshrunk * common.external_ad_control_coef > 0
    )
    common["external_bh_support"] = common.external_padj < 0.05
    deg_ids = set(deg.entrez_id) & set(common.entrez_id)
    replicated = set(
        common.loc[
            common.entrez_id.isin(deg_ids)
            & common.replicated_direction
            & common.external_bh_support,
            "entrez_id",
        ]
    )
    common["discovery_protein_coding_deg"] = common.entrez_id.isin(deg_ids)
    common["fixed_deg_externally_replicated"] = common.entrez_id.isin(replicated)
    common.to_csv(OUT / "common_gene_evidence.csv", index=False)
    return common, deg_ids, replicated


def load_reactome() -> tuple[dict[str, set[str]], dict[str, str]]:
    mapping: dict[str, set[str]] = defaultdict(set)
    with (SRC / "NCBI2Reactome_All_Levels.txt").open(encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) >= 6 and row[5] == "Homo sapiens" and row[1].startswith("R-HSA-"):
                mapping[row[1]].add(clean_id(row[0]))
    names = {}
    with (SRC / "ReactomePathways.txt").open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 3 and row[2] == "Homo sapiens":
                names[row[0]] = row[1]
    return mapping, names


def pathway_analysis(common: pd.DataFrame, deg_ids: set[str], replicated: set[str]) -> dict:
    mapping, names = load_reactome()
    all_ids = set(common.entrez_id)
    mapped_ids = set().union(*mapping.values())
    universe = all_ids & mapped_ids
    gene_symbol = dict(zip(common.entrez_id, common.symbol))
    summaries = {}
    for label, selected in (("discovery_mapped", deg_ids), ("external_replicated", replicated)):
        chosen = selected & universe
        rows = []
        for pathway_id, genes in mapping.items():
            in_universe = genes & universe
            if not 10 <= len(in_universe) <= 500:
                continue
            hits = chosen & in_universe
            p = float(hypergeom.sf(len(hits) - 1, len(universe), len(in_universe), len(chosen))) if hits else 1.0
            rows.append(
                {
                    "pathway_id": pathway_id,
                    "pathway_name": names.get(pathway_id, ""),
                    "universe_size": len(universe),
                    "pathway_universe_size": len(in_universe),
                    "selected_universe_size": len(chosen),
                    "hit_count": len(hits),
                    "pvalue": p,
                    "hit_entrez_ids": ";".join(sorted(hits, key=int)),
                    "hit_symbols": ";".join(sorted(gene_symbol.get(g, g) for g in hits)),
                }
            )
        df = pd.DataFrame(rows)
        df["padj_bh"] = multipletests(df.pvalue, method="fdr_bh")[1]
        df.sort_values(["padj_bh", "pvalue", "pathway_id"], inplace=True)
        df.to_csv(OUT / f"reactome_ora_{label}_all.csv", index=False)
        df[df.hit_count >= 3].to_csv(OUT / f"reactome_ora_{label}_presentable.csv", index=False)
        summaries[label] = {
            "input_gene_count": len(selected),
            "input_in_universe": len(chosen),
            "eligible_pathways": len(df),
            "fdr_lt_0_05_at_least_3_hits": int(((df.padj_bh < 0.05) & (df.hit_count >= 3)).sum()),
        }
    summaries["common_tested_gene_count"] = len(all_ids)
    summaries["common_reactome_universe_count"] = len(universe)
    summaries["fixed_deg_count"] = len(deg_ids)
    summaries["fixed_replicated_count"] = len(replicated)
    return summaries


def atlas_analysis(common: pd.DataFrame, deg_ids: set[str]) -> dict:
    zpath = SRC / "rna_single_nuclei_cluster_type.tsv.zip"
    with zipfile.ZipFile(zpath) as archive:
        with archive.open("rna_single_nuclei_cluster_type.tsv") as handle:
            atlas = pd.read_csv(handle, sep="\t", dtype={"Gene": str})
    atlas.rename(columns={"Gene": "ensembl_id", "Gene name": "atlas_gene_name", "Cluster type": "cluster_type", "nCPM": "ncpm"}, inplace=True)
    with zipfile.ZipFile(SRC / "rna_single_nuclei_cluster_type_cluster_types.tsv.zip") as archive:
        with archive.open("rna_single_nuclei_cluster_type_cluster_types.tsv") as handle:
            classes = pd.read_csv(handle, sep="\t")
    classes.rename(columns={"Cluster type": "cluster_type", "Cell type class": "cell_class"}, inplace=True)
    atlas = atlas.merge(classes, on="cluster_type", how="left", validate="many_to_one")
    genes = common.loc[common.entrez_id.isin(deg_ids), ["entrez_id", "symbol", "ensembl_id", "log2fc_unshrunk", "external_ad_control_coef", "external_padj", "fixed_deg_externally_replicated"]].copy()
    long = genes.merge(atlas, on="ensembl_id", how="left", validate="one_to_many")
    long["descriptively_detected_ncpm_ge_1"] = long.ncpm >= 1
    long.to_csv(OUT / "healthy_brain_cell_type_expression_long.csv", index=False)
    high = long.sort_values(["entrez_id", "ncpm"], ascending=[True, False]).drop_duplicates("entrez_id")
    high = high[["entrez_id", "symbol", "ensembl_id", "cluster_type", "cell_class", "ncpm", "fixed_deg_externally_replicated"]]
    high.rename(columns={"cluster_type": "highest_atlas_cluster_type", "cell_class": "highest_atlas_cell_class", "ncpm": "highest_atlas_ncpm"}, inplace=True)
    high.to_csv(OUT / "healthy_brain_cell_type_summary.csv", index=False)
    return {
        "fixed_degs": len(genes),
        "atlas_mapped_genes": int(high.highest_atlas_ncpm.notna().sum()),
        "atlas_unmapped_genes": int(high.highest_atlas_ncpm.isna().sum()),
        "cluster_types": int(classes.cluster_type.nunique()),
    }


def partner_atlas_analysis(common: pd.DataFrame) -> dict:
    candidates = pd.read_csv(OUT / "one_step_partner_followup_candidates.csv", dtype={"partner_entrez_id": str})
    partner_ids = set(candidates.partner_entrez_id.dropna())
    partner_genes = common.loc[common.entrez_id.isin(partner_ids), ["entrez_id", "symbol", "ensembl_id"]].copy()
    with zipfile.ZipFile(SRC / "rna_single_nuclei_cluster_type.tsv.zip") as archive:
        with archive.open("rna_single_nuclei_cluster_type.tsv") as handle:
            atlas = pd.read_csv(handle, sep="\t", dtype={"Gene": str})
    atlas.rename(columns={"Gene": "ensembl_id", "Cluster type": "cluster_type", "nCPM": "ncpm"}, inplace=True)
    long = partner_genes.merge(atlas[["ensembl_id", "cluster_type", "ncpm"]], on="ensembl_id", how="left", validate="one_to_many")
    long["descriptively_detected_ncpm_ge_1"] = long.ncpm >= 1
    long.to_csv(OUT / "healthy_brain_partner_cell_type_expression_long.csv", index=False)
    high = long.sort_values(["entrez_id", "ncpm"], ascending=[True, False]).drop_duplicates("entrez_id")
    high.rename(columns={"cluster_type": "highest_atlas_cluster_type", "ncpm": "highest_atlas_ncpm"}, inplace=True)
    high[["entrez_id", "symbol", "ensembl_id", "highest_atlas_cluster_type", "highest_atlas_ncpm"]].to_csv(
        OUT / "healthy_brain_partner_cell_type_summary.csv", index=False
    )
    return {
        "candidate_partner_genes": len(partner_genes),
        "atlas_mapped_genes": int(high.highest_atlas_ncpm.notna().sum()),
    }


ENSG = re.compile(r"ENSG\d+")


def interaction_analysis(common: pd.DataFrame, replicated: set[str]) -> dict:
    by_ensg = {str(r.ensembl_id): r for r in common.itertuples(index=False) if pd.notna(r.ensembl_id)}
    rows = []
    with (SRC / "reactome.homo_sapiens.interactions.tab-delimited.txt").open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for edge in reader:
            a = sorted(set(ENSG.findall(edge["Interactor 1 Ensembl gene id"])))
            b = sorted(set(ENSG.findall(edge["Interactor 2 Ensembl gene id"])))
            for ea in a:
                for eb in b:
                    if ea == eb or ea not in by_ensg or eb not in by_ensg:
                        continue
                    ga, gb = by_ensg[ea], by_ensg[eb]
                    for seed, partner in ((ga, gb), (gb, ga)):
                        if seed.entrez_id not in replicated:
                            continue
                        if partner.entrez_id in replicated:
                            partner_status = "fixed_deg_externally_replicated"
                        elif bool(partner.replicated_direction) and bool(partner.external_bh_support):
                            partner_status = "external_bh_and_direction;not_fixed_deg"
                        else:
                            partner_status = "insufficient_expression_support"
                        rows.append({
                            "seed_symbol": seed.symbol,
                            "seed_entrez_id": seed.entrez_id,
                            "partner_symbol": partner.symbol,
                            "partner_entrez_id": partner.entrez_id,
                            "partner_status": partner_status,
                            "partner_discovery_lfc": partner.log2fc_unshrunk,
                            "partner_external_coef": partner.external_ad_control_coef,
                            "partner_external_padj": partner.external_padj,
                            "interaction_type": edge["Interaction type"],
                            "interaction_context": edge["Interaction context"],
                            "pubmed_references": edge["Pubmed references"],
                            "source_pair_level": "Reactome pathway/complex-derived; not proof of direct binary binding",
                        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df.drop_duplicates(inplace=True)
        df.sort_values(["seed_symbol", "partner_status", "partner_symbol"], inplace=True)
    df.to_csv(OUT / "one_step_reactome_partner_evidence.csv", index=False)
    eligible = df[df.partner_status != "insufficient_expression_support"] if not df.empty else df
    eligible.to_csv(OUT / "one_step_partner_followup_candidates.csv", index=False)
    return {
        "one_step_evidence_rows": len(df),
        "one_step_unique_pairs": int(df[["seed_entrez_id", "partner_entrez_id"]].drop_duplicates().shape[0]) if not df.empty else 0,
        "expression_supported_partner_genes": int(eligible.partner_entrez_id.nunique()) if not eligible.empty else 0,
        "pair_source_limit": "Reactome-derived; direct experimental interaction unverified",
    }


if __name__ == "__main__":
    common_genes, fixed_degs, fixed_replicated = load_genes()
    summary = {"pathways": pathway_analysis(common_genes, fixed_degs, fixed_replicated)}
    summary["atlas"] = atlas_analysis(common_genes, fixed_degs)
    summary["partners"] = interaction_analysis(common_genes, fixed_replicated)
    summary["partner_atlas"] = partner_atlas_analysis(common_genes)
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
