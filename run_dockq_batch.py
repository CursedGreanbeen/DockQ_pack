#!/usr/bin/env python3
"""
Batch DockQ quality assessment for AlphaFold predictions.
Compares predicted structures against reference CIFs and outputs results to TSV.
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Optional


# Paths
FASTA_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/data/fasta-filtered")
CIFS_FILTERED_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/data/CIFs-filtered-new")
CIFS_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/data/CIFs")

# DockQ path
DOCKQ_PATH = "/home/mullagaliamova/miniconda3/envs/dockq_env/bin/DockQ"

# ANARCI report path
ANARCI_REPORT_PATH = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/results/reports/anarci_crop_report.tsv")


def load_ab_ag_mapping(report_path: Path) -> dict[str, dict[str, str]]:
    """
    Load ANARCI crop report and return chain mappings by PDB ID.
    Returns: {pdb_id: {chain: 'AB' or 'AG'}}
    """
    import csv
    from collections import defaultdict

    mappings = defaultdict(dict)

    if not report_path.exists():
        print(f"WARNING: ANARCI report not found: {report_path}")
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


def get_chain_groups(mapping: dict[str, dict[str, str]], pdb_id: str) -> tuple[list[str], list[str]]:
    """Get antibody and antigen chain lists for a PDB ID."""
    chains = mapping.get(pdb_id, {})
    ab_chains = [c for c, t in chains.items() if t == 'AB']
    ag_chains = [c for c, t in chains.items() if t == 'AG']
    return ab_chains, ag_chains


def run_command(cmd: str, capture: bool = True) -> tuple[str, int]:
    """Run a shell command and return (stdout, return_code)."""
    result = subprocess.run(
        cmd,
        shell=True,
        capture_output=capture,
        text=True
    )
    return result.stdout.strip(), result.returncode


def find_model_files(s3data_dir: Path) -> list[dict[str, Path]]:
    """
    Find all model files in s3data directory, grouped by PDB ID.
    Returns list of dicts: [{'pdb_id': '10fl', 'model_path': ..., 'ref_id': '10fl'}, ...]
    """
    models = []

    if not s3data_dir.exists():
        print(f"WARNING: s3data directory not found: {s3data_dir}")
        return models

    # First, find all top-level directories (these are the PDB IDs)
    for pdb_dir in s3data_dir.iterdir():
        if not pdb_dir.is_dir():
            continue

        pdb_id = pdb_dir.name

        for subdir in pdb_dir.iterdir():
            for model_file in subdir.rglob("*_model.cif"):
                models.append({
                    'pdb_id': pdb_id,
                    'model_path': model_file,
                    'ref_id': pdb_id  # Reference ID is the top-level directory name
                })

    return models


def find_reference(pdb_id: str, cifs_filtered_dir: Path, cifs_dir: Path) -> Optional[Path]:
    """
    Find reference structure for a PDB ID.
    First checks CIFs-filtered/, then CIFs/.
    """
    # Try CIFs-filtered first (with _cropped suffix)
    ref_path = cifs_filtered_dir / f"{pdb_id}_cropped.cif"
    if ref_path.exists():
        return ref_path

    # Try CIFs-filtered without _cropped suffix
    ref_path = cifs_filtered_dir / f"{pdb_id}.cif"
    if ref_path.exists():
        return ref_path

    # Try CIFs directory (no _cropped suffix - original files)
    ref_path = cifs_dir / f"{pdb_id}.cif"
    if ref_path.exists():
        return ref_path

    return None


def parse_dockq_output(output: str) -> tuple[dict[str, str], list[dict[str, str]]]:
    """
    Parse DockQ output to extract global metrics and per-interface metrics.
    Returns (global_results, list_of_interface_results).
    """
    global_result = {
        "dockq_score": "N/A",
    }

    interfaces = []
    current_interface = {}

    lines = output.split('\n')

    for line in lines:
        line = line.strip()

        # Parse Global DockQ score (format: "Total DockQ over X native interfaces: 0.296")
        if 'Total DockQ over' in line:
            match = re.search(r'Total DockQ over \d+ native interfaces: ([\d.]+)', line)
            if match:
                global_result["dockq_score"] = match.group(1)

        # Detect start of interface block (e.g., "Native chains: F, E")
        if line.startswith('Native chains:'):
            if current_interface:
                interfaces.append(current_interface)
            current_interface = {}
            # Extract chain info
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


def run_dockq(model_path: Path, reference_path: Path) -> Optional[tuple[dict[str, str], list[dict[str, str]]]]:
    """
    Run DockQ on a model-reference pair.
    Returns (global_results, interface_results) or None if failed.
    """
    cmd = f"{DOCKQ_PATH} {model_path} {reference_path}"
    stdout, returncode = run_command(cmd)

    if returncode != 0:
        print(f"  ERROR: DockQ failed (exit code {returncode})")
        return None

    return parse_dockq_output(stdout)


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


def aggregate_by_pdb(results: list[dict]) -> list[dict]:
    """Aggregate interface results by PDB ID, adding quality count columns."""
    pdb_data = {}

    for r in results:
        pdb_id = r.get("pdb_id", "N/A")
        if pdb_id not in pdb_data:
            pdb_data[pdb_id] = {
                "pdb_id": pdb_id,
                "model_path": r.get("model_path", "N/A"),
                "reference_path": r.get("reference_path", "N/A"),
                "dockq_scores": [],
                "iRMSD_values": [],
                "LRMSD_values": [],
                "fnat_values": [],
                "quality_counts": {"High": 0, "Medium": 0, "Acceptable": 0, "Incorrect": 0}
            }

        data = pdb_data[pdb_id]
        dockq = r.get("dockq_score", "N/A")
        if dockq != "N/A":
            data["dockq_scores"].append(float(dockq))

        irmsd = r.get("iRMSD", "N/A")
        if irmsd != "N/A":
            data["iRMSD_values"].append(float(irmsd))

        lrmsd = r.get("LRMSD", "N/A")
        if lrmsd != "N/A":
            data["LRMSD_values"].append(float(lrmsd))

        fnat = r.get("fnat", "N/A")
        if fnat != "N/A":
            data["fnat_values"].append(float(fnat))

        quality_bin = r.get("quality_bin", "N/A")
        if quality_bin in data["quality_counts"]:
            data["quality_counts"][quality_bin] += 1

    aggregated = []
    for pdb_id, data in sorted(pdb_data.items()):
        avg_dockq = f"{sum(data['dockq_scores']) / len(data['dockq_scores']):.3f}" if data["dockq_scores"] else "N/A"
        quality_bin = determine_quality_bin(avg_dockq)

        irmsd_range = f"{min(data['iRMSD_values']):.3f}-{max(data['iRMSD_values']):.3f}" if data["iRMSD_values"] else "N/A"
        lrmsd_range = f"{min(data['LRMSD_values']):.3f}-{max(data['LRMSD_values']):.3f}" if data["LRMSD_values"] else "N/A"
        fnat_range = f"{min(data['fnat_values']):.3f}-{max(data['fnat_values']):.3f}" if data["fnat_values"] else "N/A"

        aggregated.append({
            "pdb_id": pdb_id,
            "model_path": data["model_path"],
            "reference_path": data["reference_path"],
            "total_dockq": avg_dockq,
            "quality": quality_bin,
            "high_count": data["quality_counts"]["High"],
            "medium_count": data["quality_counts"]["Medium"],
            "acceptable_count": data["quality_counts"]["Acceptable"],
            "incorrect_count": data["quality_counts"]["Incorrect"],
            "iRMSD": irmsd_range,
            "LRMSD": lrmsd_range,
            "fnat": fnat_range
        })

    return aggregated


def write_tsv(results: list[dict], output_file: Path):
    """Write aggregated results to TSV file with one row per PDB."""
    headers = ["pdb_id", "Total DockQ", "Quality", "High (>=0.80)", "Medium (0.49-0.80)",
               "Acceptable (0.23-0.49)", "Incorrect (<0.23)", "iRMSD", "LRMSD", "fnat"]

    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')

        for r in results:
            row = [
                r.get("pdb_id", "N/A"),
                r.get("total_dockq", "N/A"),
                r.get("quality", "N/A"),
                str(r.get("high_count", 0)),
                str(r.get("medium_count", 0)),
                str(r.get("acceptable_count", 0)),
                str(r.get("incorrect_count", 0)),
                r.get("iRMSD", "N/A"),
                r.get("LRMSD", "N/A"),
                r.get("fnat", "N/A")
            ]
            f.write('\t'.join(row) + '\n')


def main():
    parser = argparse.ArgumentParser(description="Batch DockQ assessment for AlphaFold predictions")
    parser.add_argument(
        "--s3data-dir",
        type=Path,
        default=FASTA_DIR / "s3data",
        help=f"Directory containing AlphaFold predictions (default: {FASTA_DIR}/s3data)"
    )
    parser.add_argument(
        "--cifs-filtered-dir",
        type=Path,
        default=CIFS_FILTERED_DIR,
        help=f"Directory with filtered CIFs (default: {CIFS_FILTERED_DIR})"
    )
    parser.add_argument(
        "--cifs-dir",
        type=Path,
        default=CIFS_DIR,
        help=f"Directory with all CIFs (default: {CIFS_DIR})"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("dockq_results.tsv"),
        help="Output TSV file (default: dockq_results.tsv)"
    )
    parser.add_argument(
        "--pdb-id",
        type=str,
        help="Process only a specific PDB ID"
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="Process only the first N PDB directories"
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        help="Directory to write per-complex log files"
    )
    parser.add_argument(
        "--detailed-report",
        type=Path,
        help="Write detailed report (per-model results) to this TSV file"
    )
    parser.add_argument(
        "--anarci-report",
        type=Path,
        default=ANARCI_REPORT_PATH,
        help=f"ANARCI crop report for AB/AG chain mapping (default: {ANARCI_REPORT_PATH})"
    )
    parser.add_argument(
        "--no-chain-mapping",
        action="store_true",
        help="Skip AB/AG chain mapping even if ANARCI report is available"
    )

    args = parser.parse_args()

    # Use local variables for paths from args
    cifs_filtered_dir = args.cifs_filtered_dir
    cifs_dir = args.cifs_dir
    log_dir = args.log_dir

    # Load AB/AG chain mapping from ANARCI report
    ab_ag_mapping = {}
    if not args.no_chain_mapping and args.anarci_report.exists():
        print(f"Loading AB/AG chain mapping from: {args.anarci_report}")
        ab_ag_mapping = load_ab_ag_mapping(args.anarci_report)
        total_ab = sum(1 for m in ab_ag_mapping.values() for t in m.values() if t == 'AB')
        total_ag = sum(1 for m in ab_ag_mapping.values() for t in m.values() if t == 'AG')
        print(f"  Loaded mappings for {len(ab_ag_mapping)} complexes ({total_ab} AB chains, {total_ag} AG chains)")
    elif args.no_chain_mapping:
        print("AB/AG chain mapping disabled via --no-chain-mapping")
    else:
        print(f"WARNING: ANARCI report not found: {args.anarci_report}")

    # Create log directory if specified
    if log_dir:
        log_dir.mkdir(parents=True, exist_ok=True)

    print(f"Looking for models in: {args.s3data_dir}")
    print(f"Reference CIFs (filtered): {args.cifs_filtered_dir}")
    print(f"Reference CIFs (all): {args.cifs_dir}")
    print(f"Output file: {args.output}")
    if log_dir:
        print(f"Log directory: {log_dir}")
    if args.detailed_report:
        print(f"Detailed report: {args.detailed_report}")
    if args.limit:
        print(f"Limit: first {args.limit} PDB directories")

    # Find all model files
    if args.pdb_id:
        # Process single PDB ID - find all models under this PDB directory
        models = []
        pdb_dir = args.s3data_dir / args.pdb_id
        if pdb_dir.exists():
            for model_file in pdb_dir.rglob("*_model.cif"):
                models.append({
                    'pdb_id': args.pdb_id,
                    'model_path': model_file,
                    'ref_id': args.pdb_id
                })
        if not models:
            print(f"ERROR: No model files found for PDB ID: {args.pdb_id}")
            sys.exit(1)
    else:
        models = find_model_files(args.s3data_dir)

        # Apply limit if specified
        if args.limit:
            # Group by PDB ID first
            pdb_ids_seen = []
            limited_models = []
            for m in sorted(models, key=lambda x: x['pdb_id']):
                if m['pdb_id'] not in pdb_ids_seen:
                    pdb_ids_seen.append(m['pdb_id'])
                if len(pdb_ids_seen) <= args.limit:
                    limited_models.append(m)
            models = limited_models

    if not models:
        print("ERROR: No model files found")
        sys.exit(1)

    print(f"\nFound {len(models)} model(s)")

    # Process each model
    results = []
    success_count = 0
    fail_count = 0
    no_ref_count = 0
    pdb_logs: dict[str, list[str]] = {}  # Collect log lines per PDB

    for model_info in sorted(models, key=lambda x: (x['pdb_id'], str(x['model_path']))):
        pdb_id = model_info['pdb_id']
        model_path = model_info['model_path']
        ref_id = model_info['ref_id']

        if not model_path.exists():
            print(f"\n[{pdb_id}] Model not found: {model_path}")
            fail_count += 1
            continue

        # Find reference using ref_id (top-level directory name)
        ref_path = find_reference(ref_id, cifs_filtered_dir, cifs_dir)

        if not ref_path:
            print(f"\n[{pdb_id}] No reference found for {pdb_id}")
            no_ref_count += 1
            results.append({
                "pdb_id": pdb_id,
                "model_path": model_path,
                "reference_path": "N/A",
                "interface_id": "N/A",
                "native_chains": "N/A",
                "dockq_score": "N/A",
                "iRMSD": "N/A",
                "LRMSD": "N/A",
                "fnat": "N/A",
                "quality_bin": "No reference"
            })
            continue

        print(f"\n[{pdb_id}] Model: {model_path.name}")
        print(f"       Reference: {ref_path.name}")

        # Initialize log for this PDB if needed
        if pdb_id not in pdb_logs:
            pdb_logs[pdb_id] = []
        log_lines = pdb_logs[pdb_id]
        log_lines.append(f"Model: {model_path.name}")
        log_lines.append(f"Reference: {ref_path.name}")

        # Run DockQ
        dockq_output = run_dockq(model_path, ref_path)

        if dockq_output:
            global_result, interfaces = dockq_output
            global_dockq = global_result.get("dockq_score", "N/A")
            global_quality_bin = determine_quality_bin(global_dockq)

            print(f"       Global DockQ: {global_dockq} ({global_quality_bin})")
            print(f"       Interfaces: {len(interfaces)}")
            log_lines.append(f"Global DockQ: {global_dockq} ({global_quality_bin})")
            log_lines.append(f"Interfaces: {len(interfaces)}")

            # Create one result row per interface
            for idx, iface in enumerate(interfaces, 1):
                iface_dockq = iface.get("dockq_score", "N/A")
                iface_quality_bin = determine_quality_bin(iface_dockq)

                print(f"       Interface {idx}: DockQ={iface_dockq}, iRMSD={iface.get('iRMSD', 'N/A')}, "
                      f"LRMSD={iface.get('LRMSD', 'N/A')}, fnat={iface.get('fnat', 'N/A')} ({iface_quality_bin})")

                log_lines.append(f"  Interface {idx}: DockQ={iface_dockq}, iRMSD={iface.get('iRMSD', 'N/A')}, "
                                 f"LRMSD={iface.get('LRMSD', 'N/A')}, fnat={iface.get('fnat', 'N/A')} ({iface_quality_bin})")

                # Get AB/AG chain assignments if available
                ab_chains, ag_chains = get_chain_groups(ab_ag_mapping, pdb_id)

                results.append({
                    "pdb_id": pdb_id,
                    "model_path": model_path,
                    "reference_path": ref_path,
                    "interface_id": f"{pdb_id}_{model_path.stem}_iface{idx}",
                    "native_chains": iface.get("native_chains", "N/A"),
                    "ab_chains": ",".join(ab_chains) if ab_chains else "N/A",
                    "ag_chains": ",".join(ag_chains) if ag_chains else "N/A",
                    "dockq_score": iface_dockq,
                    "iRMSD": iface.get("iRMSD", "N/A"),
                    "LRMSD": iface.get("LRMSD", "N/A"),
                    "fnat": iface.get("fnat", "N/A"),
                    "quality_bin": iface_quality_bin
                })
            success_count += 1
        else:
            print(f"       ERROR: DockQ failed")
            log_lines.append("ERROR: DockQ failed")
            fail_count += 1

        log_lines.append("")  # Empty line between models

    # Write per-complex log files if log_dir specified
    if log_dir:
        print(f"\nWriting log files to: {log_dir}")
        for pdb_id, log_lines in pdb_logs.items():
            log_path = log_dir / f"{pdb_id}.log"
            log_path.write_text("\n".join(log_lines))
            print(f"  Written: {log_path}")

    # Write detailed report (per-model results) if requested
    if args.detailed_report:
        detailed_headers = ["pdb_id", "model_name", "reference_path", "interface_id", "native_chains",
                           "ab_chains", "ag_chains", "dockq_score", "iRMSD", "LRMSD", "fnat", "quality_bin"]
        with open(args.detailed_report, 'w') as f:
            f.write('\t'.join(detailed_headers) + '\n')
            for r in results:
                row = [
                    r.get("pdb_id", "N/A"),
                    Path(r.get("model_path", "N/A")).name,
                    str(r.get("reference_path", "N/A")),
                    r.get("interface_id", "N/A"),
                    r.get("native_chains", "N/A"),
                    r.get("ab_chains", "N/A"),
                    r.get("ag_chains", "N/A"),
                    r.get("dockq_score", "N/A"),
                    r.get("iRMSD", "N/A"),
                    r.get("LRMSD", "N/A"),
                    r.get("fnat", "N/A"),
                    r.get("quality_bin", "N/A")
                ]
                f.write('\t'.join(row) + '\n')
        print(f"\nDetailed report saved to: {args.detailed_report}")

    # Aggregate by PDB and write summary TSV
    aggregated = aggregate_by_pdb(results)
    write_tsv(aggregated, args.output)

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Total models: {len(models)}")
    print(f"Successful: {success_count}")
    print(f"No reference: {no_ref_count}")
    print(f"Failed: {fail_count}")
    print(f"\nResults saved to: {args.output}")
    print(f"Total PDB entries: {len(aggregated)}")

    # Quality distribution
    if success_count > 0:
        print(f"\nQuality distribution (per interface):")
        bins = {"High": 0, "Medium": 0, "Acceptable": 0, "Incorrect": 0}
        for r in results:
            bin_name = r.get("quality_bin", "N/A")
            if bin_name in bins:
                bins[bin_name] += 1

        for bin_name, count in bins.items():
            if count > 0:
                print(f"  {bin_name}: {count}")


if __name__ == "__main__":
    main()
