"""Chain mapping utilities for DockQ analysis."""

import csv
from collections import defaultdict
from pathlib import Path
from typing import Optional


def load_anarci_report(report_path: Path) -> dict[str, dict[str, str]]:
    """
    Load ANARCI report and return chain type mapping by PDB ID.

    Returns:
        {pdb_id: {chain: 'AB' or 'AG'}}

    Status meanings:
        - ANTIBODY -> Antibody chain (AB)
        - NOT_RECOGNIZED -> Antigen chain (AG)
    """
    mappings: dict[str, dict[str, str]] = defaultdict(dict)

    if not report_path.exists():
        return dict(mappings)

    with open(report_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            pdb_id = row['pdb_id']
            chain = row['chain']
            status = row['status']

            # Skip special status rows (but NOT NO_ANTIBODY_DOMAIN - that's antigen)
            if status == 'NO_ANTIBODY_CHAINS':
                continue

            # Strip _model suffix if present (from model ANARCI reports)
            if pdb_id.endswith('_model'):
                pdb_id = pdb_id[:-6]

            # Classify chains
            if status == 'ANTIBODY':
                mappings[pdb_id][chain] = 'AB'  # Antibody chain
            elif status in ('NOT_RECOGNIZED', 'NO_ANTIBODY_DOMAIN'):
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
