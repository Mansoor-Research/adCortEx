"""Retrieve or import the exact input snapshots listed in the source manifest.

Current-release URLs are accepted only when their bytes match the archived
SHA-256. A publisher update is an error, not permission to change the analysis.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import sys

from release_utils import sha256_file

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'config/source_manifest.json'


def prepare(check_only=False, input_dir=None, module=None):
    """Validate existing sources and, when requested, import or download missing files."""
    sources = json.loads(MANIFEST.read_text(encoding='utf-8'))['files']
    if module:
        sources = [item for item in sources if item['module'] == module]
    errors = []
    for item in sources:
        path = ROOT / item['path']
        try:
            if not path.exists():
                if check_only:
                    raise FileNotFoundError('input is missing')
                path.parent.mkdir(parents=True, exist_ok=True)
                partial = path.with_name(path.name + '.part')
                if input_dir is not None:
                    source = Path(input_dir).resolve() / item['path']
                    if not source.is_file():
                        raise FileNotFoundError('exact snapshot is absent from --input-dir')
                    shutil.copyfile(source, partial)
                else:
                    import requests
                    with requests.get(item['url'], stream=True, timeout=(20, 120)) as response:
                        response.raise_for_status()
                        with partial.open('wb') as output:
                            for chunk in response.iter_content(1024 * 1024):
                                if chunk:
                                    output.write(chunk)
                if sha256_file(partial) != item['sha256']:
                    partial.unlink()
                    raise ValueError('downloaded/imported file differs from the published analysis snapshot; obtain the exact archived file')
                partial.replace(path)
            if path.stat().st_size != item['bytes'] or sha256_file(path) != item['sha256']:
                raise ValueError('existing input does not match its recorded size and SHA-256')
            if not check_only:
                provenance = (
                    f"source_url: {item['url']}\nsha256: {item['sha256']}\n"
                    f"verified_utc: {datetime.now(timezone.utc).isoformat()}\n"
                    f"bytes: {item['bytes']}\n"
                )
                path.with_name(path.name + '.provenance.txt').write_text(provenance, encoding='utf-8')
            print('Verified:', item['path'], flush=True)
        except (OSError, ValueError) as error:
            errors.append(f"{item['path']}: {error}")
        except Exception as error:
            errors.append(f"{item['path']}: retrieval failed ({type(error).__name__}: {error})")
    if errors:
        raise RuntimeError('\n'.join(errors))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Validate inputs without downloading or writing files.')
    parser.add_argument('--input-dir', type=Path, help='Import exact snapshots from another directory with the documented data/ layout.')
    parser.add_argument('--module', choices=['discovery', 'external', 'annotations'])
    args = parser.parse_args()
    try:
        prepare(args.check, args.input_dir, args.module)
    except RuntimeError as error:
        print(error, file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
