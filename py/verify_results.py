"""Check the principal reported counts from regenerated analysis tables."""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    genes = pd.read_csv(ROOT / 'results/de/deg_protein_coding.csv')
    external = pd.read_csv(ROOT / 'results/external_validation/analysis_v2_historical_ids/gene_concordance.csv')
    overlap = pd.read_csv(ROOT / 'results/external_validation/analysis_v2_historical_ids/overlap_gene_concordance.csv')
    age = pd.read_csv(ROOT / 'results/discovery_age_sensitivity_v1/model_summary.csv')
    cells = json.loads((ROOT / 'results/cell_evidence_v1/summary.json').read_text())
    integration = json.loads((ROOT / 'results/discovery_age_sensitivity_v1/verification_and_combined_summary.json').read_text())
    categories = json.loads((ROOT / 'results/discovery_age_sensitivity_v1/integrated_category_counts.json').read_text())
    assert len(genes) == 118 and genes.entrez_id.is_unique
    assert int(external.mapped.sum()) == 112
    assert int(external.same_direction.eq(True).sum()) == 108
    assert int(external.gene_replication.sum()) == 99
    assert int(overlap.gene_replication.sum()) == 102
    assert age.set_index('model_id').fixed_supported.to_dict() == {'CC0': 118, 'CC_age': 117, 'CC_age_M3': 76}
    assert cells['pair_grid'] == 6608 and cells['finite_pair_family'] == 2214
    assert cells['significant_pairs_same'] == 389 and cells['significant_pairs_opposite'] == 181
    assert integration['combined_conditions_met'] == 40
    assert categories == {'combined_compatible_only': 29, 'combined_mixed': 11,
                          'proxy_unsupported': 34, 'proxy_supported_outside_intersection': 24,
                          'unavailable': 20}
    print('The regenerated tables reproduce the principal reported counts.')


if __name__ == '__main__':
    main()
