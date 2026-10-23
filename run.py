"""DockQ execution utilities."""

import subprocess
from pathlib import Path
from typing import Optional


DOCKQ_PATH = "/data/user/mullagaliamova/.conda/envs/dockq_env/bin/DockQ"


def run_dockq(
    model_path: Path,
    reference_path: Path,
    mappings: Optional[list[str]] = None
) -> Optional[str]:
    """
    Run DockQ on a model-reference pair.

    Args:
        model_path: Path to predicted structure
        reference_path: Path to reference structure
        mappings: Optional list of --mapping arguments (e.g., ['--mapping', 'EG:EG', '--mapping', 'FG:FG'])

    Returns:
        DockQ output as string, or None if failed
    """
    cmd = f"{DOCKQ_PATH} {model_path} {reference_path}"

    if mappings:
        cmd += ' ' + ' '.join(mappings)

    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        return None

    return result.stdout


def run_dockq_with_mappings(
    model_path: Path,
    reference_path: Path,
    ab_chains: list[str],
    ag_chains: list[str]
) -> Optional[str]:
    """
    Run DockQ without mappings. Filtering happens post-hoc based on native_chains.

    Args:
        model_path: Path to predicted structure
        reference_path: Path to reference structure
        ab_chains: List of antibody chain IDs (unused, kept for API compatibility)
        ag_chains: List of antigen chain IDs (unused, kept for API compatibility)

    Returns:
        DockQ output as string, or None if failed
    """
    # Run without --mapping, filter interfaces later
    return run_dockq(model_path, reference_path, None)
