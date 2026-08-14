"""Chain mapping utilities for DockQ analysis."""

import csv
from collections import defaultdict
from pathlib import Path
from typing import Optional


def load_anarci_report(report_path: Path) -> dict[str, dict[str, str]]:
    """
    Load ANARCI crop report and return chain type mapping by PDB ID.

    Returns:
        {pdb_id: {chain: 'AB' or 'AG'}}

    Status meanings:
        - CROPPED -> Antibody chain (AB)
        - NOT_RECOGNIZED -> Antigen chain (AG)
    """
    mappings: dict[str, dict[str, str]] = defaultdict(dict)

    if not report_path.exists():
        return dict(mappings)

    with open(report_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            file_name = row['file']
            chain = row['chain']
            status = row['status']

            # Skip special status rows
            if status in ('NO_ANTIBODY_CHAINS',):
                continue

            # Extract PDB ID from filename (e.g., "11hk.fasta" -> "11hk")
            pdb_id = file_name.replace('.fasta', '').replace('_1', '').replace('_2', '')

            # Classify chains
            if status == 'CROPPED':
                mappings[pdb_id][chain] = 'AB'  # Antibody chain
            elif status == 'NOT_RECOGNIZED':
                mappings[pdb_id][chain] = 'AG'  # Antigen chain

    return dict(mappings)


def get_chain_type_mapping(
    mapping: dict[str, dict[str, str]],
    pdb_id: str
) -> tuple[list[str], list[str]]:
    """Get antibody and antigen chain lists for a PDB ID."""
    chains = mapping.get(pdb_id, {})
    ab_chains = [c for c, t in chains.items() if t == 'AB']
    ag_chains = [c for c, t in chains.items() if t == 'AG']
    return ab_chains, ag_chains


def build_ab_ag_pairs(
    mapping: dict[str, dict[str, str]],
    pdb_id: str
) -> list[tuple[str, str]]:
    """
    Build explicit AB-AG chain pairs for DockQ --mapping flags.

    Returns list of (ab_chain, ag_chain) tuples.

    Example:
        ab_chains = ['E', 'F'], ag_chains = ['G']
        -> [('E', 'G'), ('F', 'G')]
    """
    ab_chains, ag_chains = get_chain_type_mapping(mapping, pdb_id)

    # Create all AB-AG pairs
    pairs = []
    for ab in ab_chains:
        for ag in ag_chains:
            pairs.append((ab, ag))

    return pairs


def format_dockq_mappings(pairs: list[tuple[str, str]]) -> list[str]:
    """
    Format AB-AG chain pairs as DockQ --mapping arguments.

    Syntax: --mapping MODEL_CHAINS:NATIVE_CHAINS

    Example:
        [('E', 'G'), ('F', 'G')]
        -> ['--mapping', 'EG:EG', '--mapping', 'FG:FG']
    """
    args = []
    for ab, ag in pairs:
        # Model chains: ab+ag (in order they appear in model)
        # Native chains: ab+ag (same order)
        model_chains = ab + ag
        native_chains = ab + ag
        args.extend(['--mapping', f'{model_chains}:{native_chains}'])
    return args
