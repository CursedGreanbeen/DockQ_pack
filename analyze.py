"""DockQ result analysis and reporting."""

import re
from pathlib import Path
from typing import Optional


def parse_dockq_output(output: str) -> tuple[dict[str, str], list[dict[str, str]]]:
    """
    Parse DockQ output to extract global metrics and per-interface metrics.

    Returns:
        (global_results, list_of_interface_results)
    """
    global_result = {
        "dockq_score": "N/A",
    }

    interfaces = []
    current_interface = {}

    lines = output.split('\n')

    for line in lines:
        line = line.strip()

        # Parse Global DockQ score
        if 'Total DockQ over' in line:
            match = re.search(r'Total DockQ over \d+ native interfaces: ([\d.]+)', line)
            if match:
                global_result["dockq_score"] = match.group(1)

        # Detect start of interface block
        if line.startswith('Native chains:'):
            if current_interface:
                interfaces.append(current_interface)
            current_interface = {}
            match = re.search(r'Native chains:\s*(.+)', line)
            if match:
                current_interface["native_chains"] = match.group(1)

        # Parse per-interface DockQ
        if 'DockQ:' in line and 'Total DockQ' not in line:
            match = re.search(r'DockQ:\s*([\d.]+)', line)
            if match:
                current_interface["dockq_score"] = match.group(1)

        # Parse iRMSD
        if 'iRMSD:' in line:
            match = re.search(r'iRMSD:\s*([\d.]+)', line)
            if match:
                current_interface["iRMSD"] = match.group(1)

        # Parse LRMSD
        if 'LRMSD:' in line:
            match = re.search(r'LRMSD:\s*([\d.]+)', line)
            if match:
                current_interface["LRMSD"] = match.group(1)

        # Parse fnat
        if 'fnat:' in line:
            match = re.search(r'fnat:\s*([\d.]+)', line)
            if match:
                current_interface["fnat"] = match.group(1)

    # Don't forget the last interface
    if current_interface:
        interfaces.append(current_interface)

    return global_result, interfaces


def filter_ab_ag_interfaces(
    interfaces: list[dict[str, str]],
    ab_chains: list[str],
    ag_chains: list[str]
) -> list[dict[str, str]]:
    """
    Filter interfaces to keep only antibody-antigen interfaces.

    An AB-AG interface must contain at least one AB chain AND at least one AG chain.
    Interfaces with only AB chains (e.g., E-F) or only AG chains are excluded.

    Args:
        interfaces: List of interface results from parse_dockq_output
        ab_chains: List of antibody chain IDs
        ag_chains: List of antigen chain IDs

    Returns:
        Filtered list containing only AB-AG interfaces
    """
    ab_set = set(ab_chains)
    ag_set = set(ag_chains)

    filtered = []
    for iface in interfaces:
        native_chains = iface.get("native_chains", "")
        chains = set(c.strip() for c in native_chains.replace(',', ' ').split())

        # Check if this interface has both AB and AG chains
        has_ab = bool(chains & ab_set)
        has_ag = bool(chains & ag_set)

        if has_ab and has_ag:
            filtered.append(iface)

    return filtered


def determine_quality_bin(dockq_score: str) -> str:
    """Determine quality bin from DockQ score."""
    if dockq_score == "N/A":
        return "N/A"

    try:
        score = float(dockq_score)
        if score >= 0.80:
            return "High"
        elif score >= 0.49:
            return "Medium"
        elif score >= 0.23:
            return "Acceptable"
        else:
            return "Incorrect"
    except ValueError:
        return "N/A"


def aggregate_interface_metrics(
    results: list[dict],
    group_by: str = "pdb_id"
) -> list[dict]:
    """
    Aggregate interface metrics by specified group.

    Args:
        results: List of interface result dicts
        group_by: Grouping key ("pdb_id", "ab_chains", "ab_ag_pair")

    Returns:
        Aggregated metrics per group
    """
    from collections import defaultdict

    groups: dict[str, list[dict]] = defaultdict(list)

    for r in results:
        if group_by == "pdb_id":
            key = r.get("pdb_id", "N/A")
        elif group_by == "ab_chains":
            key = r.get("ab_chains", "N/A")
        elif group_by == "ab_ag_pair":
            # Create key like "E,F:A,D,G"
            ab = r.get("ab_chains", "N/A")
            ag = r.get("ag_chains", "N/A")
            key = f"{ab}:{ag}"
        else:
            key = "all"

        groups[key].append(r)

    aggregated = []
    for key, items in sorted(groups.items()):
        dockq_scores = []
        iRMSDs = []
        fnats = []
        quality_counts = {"High": 0, "Medium": 0, "Acceptable": 0, "Incorrect": 0}

        for item in items:
            dq = item.get("dockq_score", "N/A")
            if dq != "N/A":
                dockq_scores.append(float(dq))

            ir = item.get("iRMSD", "N/A")
            if ir != "N/A":
                iRMSDs.append(float(ir))

            fn = item.get("fnat", "N/A")
            if fn != "N/A":
                fnats.append(float(fn))

            qb = item.get("quality_bin", "N/A")
            if qb in quality_counts:
                quality_counts[qb] += 1

        # Compute aggregates
        avg_dockq = f"{sum(dockq_scores) / len(dockq_scores):.3f}" if dockq_scores else "N/A"
        best_dockq = f"{max(dockq_scores):.3f}" if dockq_scores else "N/A"
        avg_iRMSD = f"{sum(iRMSDs) / len(iRMSDs):.3f}" if iRMSDs else "N/A"
        best_iRMSD = f"{min(iRMSDs):.3f}" if iRMSDs else "N/A"
        avg_fnat = f"{sum(fnats) / len(fnats):.3f}" if fnats else "N/A"
        best_fnat = f"{max(fnats):.3f}" if fnats else "N/A"

        aggregated.append({
            "group": key,
            "interface_count": len(items),
            "avg_dockq": avg_dockq,
            "best_dockq": best_dockq,
            "avg_iRMSD": avg_iRMSD,
            "best_iRMSD": best_iRMSD,
            "avg_fnat": avg_fnat,
            "best_fnat": best_fnat,
            "high_count": quality_counts["High"],
            "medium_count": quality_counts["Medium"],
            "acceptable_count": quality_counts["Acceptable"],
            "incorrect_count": quality_counts["Incorrect"],
        })

    return aggregated
