"""Published donor-pseudobulk gpath evidence for a fixed AD gene list.

Separate fetch/analyze commands preserve the outcome-join boundary.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
from urllib.parse import quote
from release_utils import git_commit

import numpy as np
import pandas as pd
import requests
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/external/cell_evidence_v1'
OUT = ROOT / 'results/cell_evidence_v1'
CONFIG = ROOT / 'config/cell_evidence_v1.json'


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_blob(content):
    return hashlib.sha1(b'blob ' + str(len(content)).encode() + b'\0' + content).hexdigest()


def get_snapshot(path, url, query=None):
    if path.exists():
        prov = json.loads(path.with_name(path.name + '.provenance.txt').read_text())
        assert sha256(path) == prov['sha256'], path
        return path
    response = requests.get(url, params=query, timeout=60)
    response.raise_for_status()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(response.content)
    prov = {'url': response.url, 'query': query, 'retrieved_utc': datetime.now(timezone.utc).isoformat(),
            'sha256': sha256(path), 'bytes': len(response.content)}
    path.with_name(path.name + '.provenance.txt').write_text(json.dumps(prov, indent=2) + '\n', encoding='utf-8')
    return path


def fetch(cfg):
    RAW.mkdir(parents=True, exist_ok=True)
    commit_file = get_snapshot(RAW / 'repository_commit.json', f"https://api.github.com/repos/{cfg['repository']}/commits/{cfg['commit']}")
    commit = json.loads(commit_file.read_text())['sha']
    assert commit == cfg['commit'], 'Published source commit does not match the fixed configuration'
    inventory = get_snapshot(RAW / 'gpath_inventory.json', f"https://api.github.com/repos/{cfg['repository']}/contents/{cfg['result_directory']}", {'ref': commit})
    files = json.loads(inventory.read_text())
    assert len(files) == cfg['expected_populations']
    assert all(x['type'] == 'file' and x['name'].endswith('_gpath.csv') for x in files)

    def one(item):
        url = f"https://raw.githubusercontent.com/{cfg['repository']}/{commit}/" + quote(cfg['result_directory'] + '/' + item['name'])
        path = get_snapshot(RAW / 'gpath' / item['name'], url)
        assert git_blob(path.read_bytes()) == item['sha'], path
        return {'filename': item['name'], 'bytes': path.stat().st_size, 'sha256': sha256(path), 'git_blob_sha1': item['sha']}

    with ThreadPoolExecutor(max_workers=6) as pool:
        rows = list(pool.map(one, files))
    manifest = {'repository': cfg['repository'], 'commit': commit, 'endpoint': cfg['endpoint'], 'files': rows}
    (RAW / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'source_commit': commit, 'files': len(rows), 'bytes': sum(r['bytes'] for r in rows)}))


def analyze(cfg):
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((RAW / 'manifest.json').read_text())
    assert len(manifest['files']) == cfg['expected_populations']
    genes = pd.read_csv(ROOT / 'results/de/deg_protein_coding.csv', dtype={'entrez_id': str})
    assert len(genes) == cfg['expected_fixed_genes'] and genes.symbol.is_unique
    comp = pd.read_csv(ROOT / 'results/de/composition_attenuation.csv', dtype={'entrez_id': str})
    comp = comp.loc[comp.model_id.eq(cfg['composition_model'])].copy()
    # Noncoding rows in the historical 156-gene sensitivity table can lack symbols.
    # Stable Entrez identity, rather than those missing labels, owns this join.
    assert comp.entrez_id.is_unique
    external = pd.read_csv(ROOT / 'results/external_validation/analysis_v2_historical_ids/gene_level_limma.csv', dtype={'entrez_id': str})
    assert external.entrez_id.is_unique
    ledger = genes.merge(external, on='entrez_id', how='left', validate='one_to_one')
    ledger = ledger.merge(comp[['entrez_id', 'adj_lfc_shrunk', 'retention', 'retention_percentile', 'adj_padj']], on='entrez_id', how='left', validate='one_to_one')
    ledger['external_mapped'] = np.isfinite(ledger.external_ad_control_coef)
    ledger['external_direction_agrees'] = ledger.external_mapped & (np.sign(ledger.external_ad_control_coef) == np.sign(ledger.log2fc_unshrunk))
    ledger['external_supported'] = ledger.external_direction_agrees & ledger.external_padj.lt(cfg['pair_bh_q'])
    ledger['m3_direction_agrees'] = np.sign(ledger.adj_lfc_shrunk) == np.sign(ledger.log2fc_unshrunk)
    ledger['m3_supported'] = ledger.m3_direction_agrees & ledger.adj_padj.lt(cfg['composition_q'])
    pairs = []
    coverage = []
    for source in manifest['files']:
        path = RAW / 'gpath' / source['filename']
        assert sha256(path) == source['sha256']
        assert git_blob(path.read_bytes()) == source['git_blob_sha1']
        frame = pd.read_csv(path)
        required = ['gene', 'cluster_id', 'logFC', 'p_val', 'p_adj.loc', 'p_adj.glb', 'coef']
        assert set(required).issubset(frame.columns), path
        assert frame.gene.notna().all() and frame.gene.is_unique, path
        assert frame.coef.eq(cfg['endpoint']).all(), path
        assert frame.cluster_id.nunique() == 1, path
        population = frame.cluster_id.iloc[0]
        for col in ['p_val', 'p_adj.loc', 'p_adj.glb']:
            values = frame[col].dropna()
            assert values.between(0, 1).all(), (path, col)
        joined = genes[['symbol', 'log2fc_unshrunk']].merge(frame[required], left_on='symbol', right_on='gene', how='left', validate='one_to_one')
        joined['source_population'] = population
        joined['source_filename'] = source['filename']
        joined['source_present'] = joined.gene.notna()
        joined['tested'] = joined.source_present & np.isfinite(joined.p_val) & np.isfinite(joined.logFC)
        joined['direction_agrees'] = joined.tested & (np.sign(joined.logFC) == np.sign(joined.log2fc_unshrunk))
        pairs.append(joined)
        coverage.append({'population': population, 'filename': source['filename'], 'source_genes': len(frame),
                         'fixed_genes_present': int(joined.source_present.sum()), 'fixed_genes_tested': int(joined.tested.sum())})
    full = pd.concat(pairs, ignore_index=True)
    assert len(full) == cfg['expected_fixed_genes'] * cfg['expected_populations']
    assert not full.duplicated(['symbol', 'source_population']).any()
    mask = full.tested
    full['module_bh_q'] = np.nan
    full['module_bonferroni_p'] = np.nan
    full.loc[mask, 'module_bh_q'] = multipletests(full.loc[mask, 'p_val'], method='fdr_bh')[1]
    full.loc[mask, 'module_bonferroni_p'] = np.minimum(1., full.loc[mask, 'p_val'] * mask.sum())
    full['pair_bh_supported'] = mask & full.module_bh_q.lt(cfg['pair_bh_q'])
    full['pair_bonferroni_supported'] = mask & full.module_bonferroni_p.lt(cfg['pair_bonferroni_p'])
    for label, support in [('bh', 'pair_bh_supported'), ('bonferroni', 'pair_bonferroni_supported')]:
        full[f'{label}_same'] = full[support] & full.direction_agrees
        full[f'{label}_opposite'] = full[support] & ~full.direction_agrees
    counts = full.groupby('symbol')[['source_present', 'tested', 'bh_same', 'bh_opposite', 'bonferroni_same', 'bonferroni_opposite']].sum().astype(int)
    ledger = ledger.merge(counts, on='symbol', validate='one_to_one')
    ledger['cell_evidence_label'] = np.select([
        ledger.tested.eq(0), ledger.bh_same.gt(0) & ledger.bh_opposite.gt(0), ledger.bh_same.gt(0), ledger.bh_opposite.gt(0)
    ], ['unavailable', 'both_directions', 'same_direction_only', 'opposite_direction_only'], default='tested_no_pair_support')
    ledger['external_m3_same_cell_support'] = ledger.external_supported & ledger.m3_supported & ledger.bh_same.gt(0)
    summary = {'fixed_genes': len(ledger), 'populations': len(coverage), 'pair_grid': len(full), 'finite_pair_family': int(mask.sum()),
               'genes_tested_any_population': int(ledger.tested.gt(0).sum()),
               'external_supported': int(ledger.external_supported.sum()), 'm3_supported': int(ledger.m3_supported.sum()),
               'genes_with_same_pair_support': int(ledger.bh_same.gt(0).sum()), 'genes_with_opposite_pair_support': int(ledger.bh_opposite.gt(0).sum()),
               'cell_labels': {str(k): int(v) for k, v in ledger.cell_evidence_label.value_counts().items()},
               'external_m3_same_cell_support': int(ledger.external_m3_same_cell_support.sum()),
               'bonferroni_same_genes': int(ledger.bonferroni_same.gt(0).sum()),
               'significant_pairs_same': int(full.bh_same.sum()), 'significant_pairs_opposite': int(full.bh_opposite.sum()),
               'source_commit': manifest['commit']}
    ledger.to_csv(OUT / 'fixed_gene_evidence_ledger.csv', index=False, lineterminator='\n')
    full.to_csv(OUT / 'fixed_gene_population_evidence.csv', index=False, lineterminator='\n')
    pd.DataFrame(coverage).to_csv(OUT / 'population_mapping_coverage.csv', index=False, lineterminator='\n')
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    inputs = ['config/cell_evidence_v1.json', 'results/de/deg_protein_coding.csv',
              'results/de/composition_attenuation.csv', 'results/external_validation/analysis_v2_historical_ids/gene_level_limma.csv',
              'data/external/cell_evidence_v1/manifest.json', 'py/21_cell_evidence.py']
    run = {'run_utc': datetime.now(timezone.utc).isoformat(), 'git_commit': git_commit(ROOT),
           'python': platform.python_version(), 'pandas': pd.__version__, 'numpy': np.__version__,
           'inputs_sha256': {name: sha256(ROOT / name) for name in inputs}}
    (OUT / 'run_info.json').write_text(json.dumps(run, indent=2) + '\n', encoding='utf-8')
    session = ROOT / 'env/session_info' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    session.mkdir(parents=True, exist_ok=True)
    (session / 'cell_evidence_run.json').write_text(json.dumps(run, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['fetch', 'analyze'])
    args = parser.parse_args()
    configuration = json.loads(CONFIG.read_text())
    {'fetch': fetch, 'analyze': analyze}[args.mode](configuration)
