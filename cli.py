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
from dockq.ranking import filter_models_by_ranking
from dockq.run import run_dockq, run_dockq_with_mappings
from dockq.analyze import (
    parse_dockq_output,
    filter_ab_ag_interfaces,
    determine_quality_bin,
    aggregate_interface_metrics,
    is_ab_ag_interface,
)


# Default paths
FASTA_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/data/fasta-filtered")
CIFS_FILTERED_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/data/CIFs-filtered-new")
CIFS_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/data/CIFs")
ANARCI_REPORT_PATH = Path("/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/reports/anarci_report.tsv")


def find_model_files(s3data_dir: Path, pdb_id: str = None) -> list[dict[str, Path]]:
    """
    Find all model files in s3data directory.

    For multiseed layouts (<pdb>/seed-<seed>_sample-<k>/), only the seed
    models are returned; the top-level <pdb>_model.cif (the pipeline's
    overall best) is skipped. For single-model layouts (no seed dirs),
    the top-level model is returned.

    Args:
        s3data_dir: Directory containing AlphaFold predictions
        pdb_id: If set, only process this PDB ID
    """
    models = []

    if not s3data_dir.exists():
        print(f"WARNING: s3data directory not found: {s3data_dir}")
        return models

    pdb_dirs = [d for d in s3data_dir.iterdir() if d.is_dir()]
    if pdb_id is not None:
        pdb_dirs = [d for d in pdb_dirs if d.name == pdb_id]

    for pdb_dir in pdb_dirs:
        current_pdb = pdb_dir.name

        # Find seed subdirectories
        seed_dirs = [d for d in pdb_dir.iterdir() if d.is_dir() and d.name.startswith('seed-')]

        if seed_dirs:
            # Multiseed: one model per seed-<seed>_sample-<k> directory
            for seed_dir in sorted(seed_dirs):
                for model_file in seed_dir.rglob("*_model.cif"):
                    models.append({
                        'pdb_id': current_pdb,
                        'model_path': model_file,
                        'ref_id': current_pdb,
                        'seed': seed_dir.name
                    })
        else:
            # No seed subdirectories - single model, find it directly
            for model_file in pdb_dir.rglob("*_model.cif"):
                models.append({
                    'pdb_id': current_pdb,
                    'model_path': model_file,
                    'ref_id': current_pdb
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


def write_header(output_file: Path):
    """Write TSV header if file doesn't exist."""
    headers = [
        "pdb_id", "model_name", "reference_path", "interface_id",
        "native_chains", "ab_chains", "ag_chains",
        "dockq_score", "iRMSD", "LRMSD", "fnat", "quality_bin"
    ]
    if not output_file.exists() or output_file.stat().st_size == 0:
        with open(output_file, 'w') as f:
            f.write('\t'.join(headers) + '\n')


def append_results(results: list[dict], output_file: Path):
    """Append results to TSV file."""
    with open(output_file, 'a') as f:
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


def write_interface_report(results: list[dict], output_file: Path):
    """
    Write interface report with AB-AG interface metrics only.

    One row per interface, filtered to include only antibody-antigen interfaces.
    """
    write_header(output_file)
    append_results(results, output_file)
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
        default=Path("dockq_detailed_report.tsv"),
        help="Output TSV file with all interfaces (default: dockq_detailed_report.tsv)"
    )
    parser.add_argument(
        "--interface-report",
        type=Path,
        help="Write AB-AG interfaces only to this file (filtered report)"
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
        "--top-per-seed",
        type=int,
        default=None,
        help="Only process the top N models per seed by ranking_score "
             "(needs <pdb>_ranking_scores.csv next to the seed dirs)"
    )
    parser.add_argument(
        "--top-per-seed",
        type=int,
        default=None,
        help="Only process the top N models per seed by ranking_score "
             "(needs <pdb>_ranking_scores.csv next to the seed dirs)"
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
        help=f"ANARCI report for AB/AG chain mapping (default: {ANARCI_REPORT_PATH})"
    )
    parser.add_argument(
        "--no-chain-mapping",
        action="store_true",
        help="Skip AB/AG chain mapping even if ANARCI report is available"
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
    print(f"Output file (all interfaces): {args.output}")
    if args.interface_report:
        print(f"Interface report (AB-AG only): {args.interface_report}")
    if args.summary:
        print(f"Summary report: {args.summary}")

    # Find all model files
    if args.pdb_id:
        models = find_model_files(args.s3data_dir, pdb_id=args.pdb_id)
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

    # Pre-filter to top N models per seed (needs <pdb>_ranking_scores.csv)
    if args.top_per_seed is not None:
        models = filter_models_by_ranking(models, args.top_per_seed)
        print(f"After ranking filter (top {args.top_per_seed}/seed): {len(models)} model(s)")

    if not models:
        print("ERROR: No model files found")
        sys.exit(1)

    print(f"\nFound {len(models)} model(s)")

    # Initialize output files with headers
    write_header(args.output)
    if args.interface_report:
        write_header(args.interface_report)

    # Process each model
    all_results = []  # Keep for summary report
    success_count = 0
    fail_count = 0
    no_ref_count = 0

    for model_info in sorted(models, key=lambda x: (x['pdb_id'], str(x['model_path']))):
        pdb_id = model_info['pdb_id']
        model_path = model_info['model_path']
        ref_id = model_info['ref_id']

        # Reset per-model results list
        results = []

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

        # Get chain mapping for filtering AB-AG interfaces
        ab_chains, ag_chains = get_chain_type_mapping(ab_ag_mapping, pdb_id)

        # Run DockQ (no --mapping flags, filter by native_chains post-hoc)
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

        # Create results for each interface (all interfaces, not filtered)
        if ab_chains and ag_chains:
            ab_ag_count = sum(1 for iface in interfaces if is_ab_ag_interface(iface.get("native_chains", ""), ab_chains, ag_chains))
            if ab_ag_count < len(interfaces):
                print(f"       AB-AG interfaces: {ab_ag_count} (excluded {len(interfaces) - ab_ag_count} non-AB-AG)")

        for idx, iface in enumerate(interfaces, 1):
            iface_dockq = iface.get("dockq_score", "N/A")
            quality_bin = determine_quality_bin(iface_dockq)

            print(f"       Interface {idx}: {iface.get('native_chains', 'N/A')} "
                  f"DockQ={iface_dockq} ({quality_bin})")

            native_chains = iface.get("native_chains", "N/A")
            is_ab_ag = is_ab_ag_interface(native_chains, ab_chains, ag_chains) if ab_chains and ag_chains else True

            result = {
                "pdb_id": pdb_id,
                "model_path": model_path,
                "reference_path": ref_path,
                "interface_id": f"{pdb_id}_{model_path.stem}_iface{idx}",
                "native_chains": native_chains,
                "ab_chains": ",".join(ab_chains) if ab_chains else "N/A",
                "ag_chains": ",".join(ag_chains) if ag_chains else "N/A",
                "dockq_score": iface_dockq,
                "iRMSD": iface.get("iRMSD", "N/A"),
                "LRMSD": iface.get("LRMSD", "N/A"),
                "fnat": iface.get("fnat", "N/A"),
                "quality_bin": quality_bin,
                "is_ab_ag": is_ab_ag,
            }
            results.append(result)
            all_results.append(result)

        if dockq_output is not None:
            # Write results incrementally after each model
            append_results(results, args.output)
            if args.interface_report:
                ab_ag_results = [r for r in results if r.get("is_ab_ag", False)]
                append_results(ab_ag_results, args.interface_report)

        success_count += 1

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Total models: {len(models)}")
    print(f"Successful: {success_count}")
    print(f"No reference: {no_ref_count}")
    print(f"Failed: {fail_count}")
    print(f"Total interfaces: {len(all_results)}")


if __name__ == "__main__":
    main()
