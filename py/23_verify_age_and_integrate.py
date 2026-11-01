"""Check fixed-list age sensitivity and integrate its descriptive evidence ledger."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, norm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/discovery_age_sensitivity_v1'
cfg = json.loads((ROOT / 'config/discovery_age_sensitivity_v1.json').read_text())
hashes = json.loads((OUT / 'input_hashes.json').read_text())
for name, expected in hashes.items():
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
fixed = pd.read_csv(OUT / 'fixed_gene_age_sensitivity.csv', dtype={'entrez_id': str})
assert len(fixed) == cfg['expected_fixed_genes'] * len(cfg['models'])
ledger = pd.read_csv(ROOT / 'results/cell_evidence_v1/fixed_gene_evidence_ledger.csv', dtype={'entrez_id': str})
checks = []
for model in cfg['models']:
    all_genes = pd.read_csv(OUT / f'{model}_all_genes.csv', dtype={'entrez_id': str})
    diagnostics = pd.read_csv(OUT / f'{model}_fit_diagnostics.csv', dtype={0: str})
    assert len(all_genes) == cfg['expected_genes'] and all_genes.entrez_id.is_unique
    chosen = fixed.loc[fixed.model_id.eq(model)].copy()
    lookup = all_genes.set_index('entrez_id').loc[chosen.entrez_id]
    assert np.allclose(lookup.log2FoldChange, chosen.log2FoldChange)
    assert np.allclose(lookup.padj, chosen.padj_sensitivity, equal_nan=True)
    valid = all_genes.padj.notna()
    assert np.allclose(false_discovery_control(all_genes.loc[valid, 'pvalue'].to_numpy(), method='bh'), all_genes.loc[valid, 'padj'], atol=1e-12)
    finite = all_genes.pvalue.notna()
    assert np.allclose(2 * norm.sf(np.abs(all_genes.loc[finite, 'stat'])), all_genes.loc[finite, 'pvalue'], atol=1e-12)
    direction = np.sign(chosen.log2FoldChange) == np.sign(chosen.log2fc_unshrunk)
    assert np.array_equal(direction, chosen.direction_agrees)
    assert np.array_equal(direction & chosen.padj_sensitivity.lt(cfg['q_threshold']), chosen.supported)
    failed_ids = diagnostics.loc[diagnostics.betaConv.eq(False)].iloc[:, 0].astype(str)
    assert not chosen.entrez_id.isin(failed_ids).any()
    # Disclose background nonconvergence; no fixed-list gene uses such a fit.
    checks.append({'model': model, 'fixed_supported': int(chosen.supported.sum()), 'fixed_direction_agrees': int(direction.sum()),
                   'background_nonconverged_ids': failed_ids.tolist(), 'fixed_nonconverged': 0,
                   'q_unavailable_all_genes': int(all_genes.padj.isna().sum())})
    part = chosen[['entrez_id', 'log2FoldChange', 'padj_sensitivity', 'direction_agrees', 'supported']].copy()
    part.columns = ['entrez_id'] + [model + '_' + x for x in ['lfc', 'q', 'direction_agrees', 'supported']]
    ledger = ledger.merge(part, on='entrez_id', validate='one_to_one')
ledger['combined_conditions_met'] = (ledger.external_supported & ledger.m3_supported & ledger.CC_age_supported &
                                     ledger.CC_age_M3_supported & ledger.bh_same.gt(0))
ledger.to_csv(OUT / 'combined_fixed_gene_evidence.csv', index=False, lineterminator='\n')
combined = ledger.loc[ledger.combined_conditions_met]
result = {'models_verified': checks, 'fixed_genes': len(ledger), 'combined_conditions_met': len(combined),
          'combined_same_direction_only': int(combined.cell_evidence_label.eq('same_direction_only').sum()),
          'combined_both_directions': int(combined.cell_evidence_label.eq('both_directions').sum()),
          'external_plus_age_M3_supported': int((ledger.external_supported & ledger.CC_age_M3_supported).sum()),
          'statement': 'Descriptive intersection of correlated evidence conditions; not a validated target score or gene-level FDR-controlled selection.'}
(OUT / 'verification_and_combined_summary.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))

# Assign the five descriptive manuscript categories without applying a new selection rule.
category = np.select([
    ledger.tested.eq(0),
    ledger.combined_conditions_met & ledger.cell_evidence_label.eq('same_direction_only'),
    ledger.combined_conditions_met & ledger.cell_evidence_label.eq('both_directions'),
    ~ledger.CC_age_M3_supported,
], ['unavailable', 'combined_compatible_only', 'combined_mixed', 'proxy_unsupported'],
    default='proxy_supported_outside_intersection')
categories = ledger[['entrez_id', 'symbol']].copy()
categories['evidence_category'] = category
categories.to_csv(OUT / 'integrated_evidence_categories.csv', index=False)
category_counts = categories.evidence_category.value_counts().to_dict()
(OUT / 'integrated_category_counts.json').write_text(json.dumps(category_counts, indent=2) + '\n', encoding='utf-8')
