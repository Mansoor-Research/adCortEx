"""Integrity checks for frozen sources and the biological-evidence extension."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
SOURCES = ROOT / "data/external/biological_evidence"
OUT = ROOT / "results/biological_evidence"

for path in SOURCES.rglob("*"):
    if not path.is_file():
        continue
    if path.name.endswith(".provenance.txt"):
        continue
    provenance = path.with_name(path.name + ".provenance.txt").read_text(encoding="utf-8")
    expected = next(line.split(": ", 1)[1] for line in provenance.splitlines() if line.startswith("sha256: "))
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    assert actual == expected, f"Source SHA mismatch: {path.name}"

for label in ("discovery_mapped", "external_replicated"):
    table = pd.read_csv(OUT / f"reactome_ora_{label}_all.csv")
    assert len(table) == 1641 and table.pathway_id.is_unique
    recomputed = hypergeom.sf(
        table.hit_count.to_numpy() - 1,
        table.universe_size.to_numpy(),
        table.pathway_universe_size.to_numpy(),
        table.selected_universe_size.to_numpy(),
    )
    assert max(abs(recomputed - table.pvalue.to_numpy())) < 1e-12
    assert max(abs(multipletests(table.pvalue, method="fdr_bh")[1] - table.padj_bh)) < 1e-12

atlas = pd.read_csv(OUT / "healthy_brain_cell_type_summary.csv")
assert len(atlas) == 112 and atlas.entrez_id.is_unique
partners = pd.read_csv(OUT / "one_step_partner_followup_candidates.csv")
assert set(partners.partner_status) == {
    "fixed_deg_externally_replicated", "external_bh_and_direction;not_fixed_deg"
}
print("Biological-evidence source, ORA, atlas and partner checks passed")
