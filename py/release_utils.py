"""Small provenance utilities shared by the public analysis scripts."""
from pathlib import Path
import hashlib
import subprocess


def sha256_file(path: Path) -> str:
    """Hash a file in bounded-memory blocks."""
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def git_commit(root: Path) -> str | None:
    """Return a revision only when this directory is itself the Git repository root.

    Downloaded ZIP archives are valid analysis inputs. A parent repository must
    not be mistaken for the public code repository.
    """
    try:
        top = subprocess.check_output(
            ['git', 'rev-parse', '--show-toplevel'], cwd=root,
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
        if Path(top).resolve() != root.resolve():
            return None
        return subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=root,
            text=True, stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None
