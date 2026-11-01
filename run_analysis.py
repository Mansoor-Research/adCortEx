"""Run the analysis-only cortical Alzheimer's disease workflow.

Every subprocess uses project-relative inputs and writes to ignored data/,
results/ or env/ directories. No figure or manuscript generation is included.
"""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'py'))
from fetch_sources import prepare

STAGES = [
    ('discovery_inputs', 'R', 'R/01_data_acquisition.R', []),
    ('discovery_qc', 'R', 'R/02_qc_covariates.R', []),
    ('discovery_de', 'R', 'R/03_differential_expression.R', []),
    ('fixed_genes', 'R', 'R/04_protein_coding_filter.R', []),
    ('composition', 'R', 'R/03b_cell_composition.R', []),
    ('external_metadata', 'Python', 'py/11_prepare_gse33000_metadata.py', []),
    ('external', 'R', 'R/11_external_validation.R', []),
    ('external_overlap', 'R', 'R/11_external_overlap_sensitivity.R', []),
    ('annotations', 'Python', 'py/12_biological_evidence.py', []),
    ('verify_annotations', 'Python', 'py/12_verify_biological_evidence.py', []),
    ('cell_fetch', 'Python', 'py/21_cell_evidence.py', ['fetch']),
    ('cell_evidence', 'Python', 'py/21_cell_evidence.py', ['analyze']),
    ('age_sensitivity', 'R', 'R/12_discovery_age_sensitivity.R', []),
    ('integration', 'Python', 'py/23_verify_age_and_integrate.py', []),
    ('verify_results', 'Python', 'py/verify_results.py', []),
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch', action='store_true', help='Retrieve missing source snapshots before running.')
    parser.add_argument('--input-dir', type=Path, help='Import missing snapshots from a directory with the documented data/ layout.')
    parser.add_argument('--check', action='store_true', help='Check dependencies and source inputs without fitting models.')
    parser.add_argument('--rscript', default='Rscript', help='Rscript executable or its full path.')
    parser.add_argument('--steps', nargs='+', choices=[item[0] for item in STAGES], help='Run selected stages in dependency order; prior outputs must already exist.')
    args = parser.parse_args()
    rscript = shutil.which(args.rscript)
    if rscript is None:
        parser.error('Rscript was not found. Add R to PATH or provide --rscript.')
    try:
        import numpy, pandas, requests, scipy, statsmodels
        subprocess.run([rscript, str(ROOT / 'R/00_check_environment.R')], cwd=ROOT, check=True)
        if args.check:
            prepare(check_only=True)
            print('Dependencies and fixed source inputs are available. No models were fitted.')
            return 0
        prepare(check_only=not (args.fetch or args.input_dir), input_dir=args.input_dir)
        for label, language, script, extra in STAGES:
            if args.steps and label not in args.steps:
                continue
            print(f'Running {label}', flush=True)
            executable = rscript if language == 'R' else sys.executable
            subprocess.run([executable, str(ROOT / script), *extra], cwd=ROOT, check=True)
    except (ImportError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f'Analysis stopped: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
