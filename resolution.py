"""Add resolution data to DockQ results from SAbDab summary."""

import argparse
import sys
from pathlib import Path
from typing import Optional


def load_sabdab_resolutions(sabdab_path: Path) -> dict[str, float]:
    """
    Load resolution data from SAbDab summary file.

    Args:
        sabdab_path: Path to SAbDab TSV file (e.g., sabdab_2024_hr_full_final.tsv)

    Returns:
        Dictionary mapping PDB ID → resolution (in Angstroms)
    """
    resolutions = {}

    with open(sabdab_path, 'r') as f:
        header = f.readline().strip().split('\t')

        pdb_idx = header.index('pdb')
        res_idx = header.index('resolution')

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > max(pdb_idx, res_idx):
                pdb_id = parts[pdb_idx]
                try:
                    resolution = float(parts[res_idx])
                    resolutions[pdb_id] = resolution
                except ValueError:
                    continue

    return resolutions


def add_resolution_to_dockq_results(
    dockq_path: Path,
    resolutions: dict[str, float],
    output_path: Path
) -> int:
    """
    Add resolution column to DockQ results file.

    Args:
        dockq_path: Path to input DockQ TSV file
        resolutions: Dictionary mapping PDB ID → resolution
        output_path: Path to output TSV file

    Returns:
        Number of complexes with resolution data added
    """
    with open(dockq_path, 'r') as f:
        lines = f.readlines()

    if not lines:
        print("ERROR: Empty input file")
        return 0

    header = lines[0].strip().split('\t')
    pdb_idx = header.index('pdb_id')

    # Add resolution column
    new_header = header + ['resolution']
    matched_count = 0

    with open(output_path, 'w') as f:
        f.write('\t'.join(new_header) + '\n')

        for line in lines[1:]:
            parts = line.strip().split('\t')
            pdb_id = parts[pdb_idx] if len(parts) > pdb_idx else None

            if pdb_id and pdb_id in resolutions:
                resolution = f"{resolutions[pdb_id]:.2f}"
                matched_count += 1
            else:
                resolution = "N/A"

            f.write('\t'.join(parts) + '\t' + resolution + '\n')

    return matched_count


def main():
    parser = argparse.ArgumentParser(
        description="Add resolution data to DockQ results from SAbDab"
    )
    parser.add_argument(
        "--sabdab",
        type=Path,
        default=Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/FOLDING_DATA/sabdab_summary_filtering/sabdab_2024_hr_full_final.tsv"),
        help="Path to SAbDab summary file (default: sabdab_2024_hr_full_final.tsv)"
    )
    parser.add_argument(
        "--dockq",
        type=Path,
        required=True,
        help="Path to DockQ results TSV file"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        help="Output TSV file (default: <input>_with_resolution.tsv)"
    )

    args = parser.parse_args()

    # Load resolution data
    print(f"Loading resolutions from: {args.sabdab}")
    resolutions = load_sabdab_resolutions(args.sabdab)
    print(f"  Loaded {len(resolutions)} resolution entries")

    # Determine output path
    if args.output:
        output_path = args.output
    else:
        # Default to results/reports directory
        RESULTS_REPORTS_DIR = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/results/reports")
        RESULTS_REPORTS_DIR.mkdir(parents=True, exist_ok=True)
        output_path = RESULTS_REPORTS_DIR / f"{args.dockq.stem}_with_resolution.tsv"

    # Add resolution to results
    print(f"Processing: {args.dockq}")
    matched = add_resolution_to_dockq_results(args.dockq, resolutions, output_path)

    print(f"Output saved to: {output_path}")
    print(f"Complexes with resolution: {matched}")
    print(f"Complexes without resolution: {len(open(args.dockq).readlines()) - 1 - matched}")


if __name__ == "__main__":
    main()
