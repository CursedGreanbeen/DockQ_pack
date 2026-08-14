#!/usr/bin/env python3
"""
Batch DockQ quality assessment for AlphaFold predictions.
Analyzes antibody-antigen interfaces with explicit chain mappings.
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from dockq.mapping import (
    load_anarci_report,
    get_chain_type_mapping,
    build_ab_ag_pairs,
    format_dockq_mappings,
)
from dockq.run import run_dockq, run_dockq_with_mappings
from dockq.analyze import (
    parse_dockq_output,
    filter_ab_ag_interfaces,
    determine_quality_bin,
    aggregate_interface_metrics,
)


# Default paths
FASTA_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/data/fasta-filtered")
CIFS_FILTERED_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/data/CIFs-filtered-new")
CIFS_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/data/CIFs")
ANARCI_REPORT_PATH = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/results/reports/anarci_crop_report.tsv")


def find_model_files(s3data_dir: Path) -> list[dict[str, Path]]:
    """Find all model files in s3data directory."""
    models = []

    if not s3data_dir.exists():
        print(f"WARNING: s3data directory not found: {s3data_dir}")
        return models

    for pdb_dir in s3data_dir.iterdir():
        if not pdb_dir.is_dir():
            continue

        pdb_id = pdb_dir.name

        for model_file in pdb_dir.rglob("*_model.cif"):
            models.append({
                'pdb_id': pdb_id,
                'model_path': model_file,
                'ref_id': pdb_id
            })

    return models


def find_reference(pdb_id: str, cifs_filtered_dir: Path, cifs_dir: Path) -> Optional[Path]:
    """Find reference structure for a PDB ID."""
    for base_dir in [cifs_filtered_dir, cifs_dir]:
        for suffix in ["_cropped.cif", ".cif"]:
            ref_path = base_dir / f"{pdb_id}{suffix}"
            if ref_path.exists():
                return ref_path

    return None


def write_interface_report(results: list[dict], output_file: Path):
    """
    Write interface report with AB-AG interface metrics only.

    One row per interface, filtered to include only antibody-antigen interfaces.
    """
    headers = [
        "pdb_id", "model_name", "reference_path", "interface_id",
        "native_chains", "ab_chains", "ag_chains",
        "dockq_score", "iRMSD", "LRMSD", "fnat", "quality_bin"
    ]

    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')

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
                r.get("quality_bin", "N/A"),
            ]
            f.write('\t'.join(row) + '\n')

    print(f"Interface report saved to: {output_file}")


def write_summary_report(aggregated: list[dict], output_file: Path, group_by: str):
    """Write summary report with aggregated metrics."""
    headers = [
        group_by, "interface_count",
        "best_dockq", "avg_dockq",
        "best_iRMSD", "avg_iRMSD",
        "best_fnat", "avg_fnat",
        "high_count", "medium_count", "incorrect_count"
    ]

    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')

        for r in aggregated:
            row = [
                r.get("group", "N/A"),
                str(r.get("interface_count", 0)),
                r.get("best_dockq", "N/A"),
                r.get("avg_dockq", "N/A"),
                r.get("best_iRMSD", "N/A"),
                r.get("avg_iRMSD", "N/A"),
                r.get("best_fnat", "N/A"),
                r.get("avg_fnat", "N/A"),
                str(r.get("high_count", 0)),
                str(r.get("medium_count", 0)),
                str(r.get("incorrect_count", 0)),
            ]
            f.write('\t'.join(row) + '\n')

    print(f"Summary report saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Batch DockQ assessment for antibody-antigen complexes"
    )
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
        default=Path("dockq_interface_results.tsv"),
        help="Output TSV file for interface report (default: dockq_interface_results.tsv)"
    )
    parser.add_argument(
        "--summary",
        type=Path,
        help="Write summary report to this file"
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
    parser.add_argument(
        "--use-mappings",
        action="store_true",
        help="Pass explicit --mapping flags to DockQ for AB-AG pairs only"
    )
    parser.add_argument(
        "--group-by",
        type=str,
        choices=["pdb_id", "ab_chains", "ab_ag_pair"],
        default="pdb_id",
        help="Grouping key for summary report (default: pdb_id)"
    )

    args = parser.parse_args()

    # Load AB/AG chain mapping from ANARCI report
    ab_ag_mapping = {}
    if not args.no_chain_mapping and args.anarci_report.exists():
        print(f"Loading AB/AG chain mapping from: {args.anarci_report}")
        ab_ag_mapping = load_anarci_report(args.anarci_report)
        total_ab = sum(1 for m in ab_ag_mapping.values() for t in m.values() if t == 'AB')
        total_ag = sum(1 for m in ab_ag_mapping.values() for t in m.values() if t == 'AG')
        print(f"  Loaded mappings for {len(ab_ag_mapping)} complexes ({total_ab} AB chains, {total_ag} AG chains)")
    elif args.no_chain_mapping:
        print("AB/AG chain mapping disabled via --no-chain-mapping")
    else:
        print(f"WARNING: ANARCI report not found: {args.anarci_report}")

    # Create log directory if specified
    if args.log_dir:
        args.log_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nLooking for models in: {args.s3data_dir}")
    print(f"Reference CIFs (filtered): {args.cifs_filtered_dir}")
    print(f"Reference CIFs (all): {args.cifs_dir}")
    print(f"Output file: {args.output}")
    if args.summary:
        print(f"Summary report: {args.summary}")
    if args.use_mappings:
        print("Using explicit --mapping flags for AB-AG pairs only")

    # Find all model files
    if args.pdb_id:
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

        if args.limit:
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

    for model_info in sorted(models, key=lambda x: (x['pdb_id'], str(x['model_path']))):
        pdb_id = model_info['pdb_id']
        model_path = model_info['model_path']
        ref_id = model_info['ref_id']

        if not model_path.exists():
            print(f"\n[{pdb_id}] Model not found: {model_path}")
            fail_count += 1
            continue

        # Find reference
        ref_path = find_reference(ref_id, args.cifs_filtered_dir, args.cifs_dir)

        if not ref_path:
            print(f"\n[{pdb_id}] No reference found for {pdb_id}")
            no_ref_count += 1
            continue

        print(f"\n[{pdb_id}] Model: {model_path.name}")
        print(f"       Reference: {ref_path.name}")

        # Get chain mapping
        ab_chains, ag_chains = get_chain_type_mapping(ab_ag_mapping, pdb_id)

        if args.use_mappings and ab_chains and ag_chains:
            # Run DockQ with explicit AB-AG mappings
            pairs = build_ab_ag_pairs(ab_ag_mapping, pdb_id)
            mapping_args = format_dockq_mappings(pairs)
            print(f"       Mappings: {pairs}")
            dockq_output = run_dockq(model_path, ref_path, mapping_args)
        else:
            # Run DockQ without explicit mappings
            dockq_output = run_dockq(model_path, ref_path)

        if dockq_output is None:
            print(f"       ERROR: DockQ failed")
            fail_count += 1
            continue

        # Parse output
        global_result, interfaces = parse_dockq_output(dockq_output)
        global_dockq = global_result.get("dockq_score", "N/A")
        print(f"       Global DockQ: {global_dockq}")
        print(f"       Interfaces found: {len(interfaces)}")

        # Filter to AB-AG interfaces only
        if ab_chains and ag_chains:
            filtered = filter_ab_ag_interfaces(interfaces, ab_chains, ag_chains)
            if len(filtered) < len(interfaces):
                print(f"       AB-AG interfaces: {len(filtered)} (excluded {len(interfaces) - len(filtered)} non-AB-AG)")
            interfaces = filtered

        # Create results for each interface
        for idx, iface in enumerate(interfaces, 1):
            iface_dockq = iface.get("dockq_score", "N/A")
            quality_bin = determine_quality_bin(iface_dockq)

            print(f"       Interface {idx}: {iface.get('native_chains', 'N/A')} "
                  f"DockQ={iface_dockq} ({quality_bin})")

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
                "quality_bin": quality_bin,
            })

        success_count += 1

    # Write interface report
    write_interface_report(results, args.output)

    # Write summary report if requested
    if args.summary:
        aggregated = aggregate_interface_metrics(results, group_by=args.group_by)
        write_summary_report(aggregated, args.summary, args.group_by)

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Total models: {len(models)}")
    print(f"Successful: {success_count}")
    print(f"No reference: {no_ref_count}")
    print(f"Failed: {fail_count}")
    print(f"Total interfaces: {len(results)}")


if __name__ == "__main__":
    main()
