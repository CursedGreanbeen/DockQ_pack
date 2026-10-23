#!/usr/bin/env python3
"""Recreate an AB-AG interface-only report from a full per-interface detailed report.

DockQ names interface chains (native_chains) with the REFERENCE/auth chain IDs, while the
ab_chains/ag_chains columns of the detailed report are in the MODEL chain namespace.
Matching the two namespaces is what truncated the old interface reports. This script
filters against anarci_report_ref.tsv (reference namespace) instead, and rewrites the
ab_chains/ag_chains columns into the reference namespace so each output row is
self-consistent with its native_chains.

Usage: recreate_interface_report_ref.py <detailed.tsv> <interface.tsv> [anarci_ref.tsv]
"""
import csv, sys
from collections import defaultdict
from pathlib import Path

STATUS_AB = {"ANTIBODY"}
STATUS_AG = {"NO_ANTIBODY_DOMAIN", "NOT_RECOGNIZED"}

DEFAULT_REF = "/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/results/reports/anarci_report_ref.tsv"

detailed = Path(sys.argv[1])
dst = Path(sys.argv[2])
anarci_ref = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(DEFAULT_REF)

def chain_set(s):
    return {c.strip() for c in s.replace(",", " ").split() if c.strip() and c.strip() != "N/A"}

ab_map, ag_map = defaultdict(set), defaultdict(set)
with open(anarci_ref, newline="") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["status"] in STATUS_AB:
            ab_map[row["pdb_id"]].add(row["chain"])
        elif row["status"] in STATUS_AG:
            ag_map[row["pdb_id"]].add(row["chain"])

kept = total = unmapped = 0
with open(detailed, newline="") as fi, open(dst, "w", newline="") as fo:
    reader = csv.DictReader(fi, delimiter="\t")
    writer = csv.DictWriter(fo, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in reader:
        total += 1
        ab, ag = ab_map.get(row["pdb_id"], set()), ag_map.get(row["pdb_id"], set())
        if not (ab and ag):
            # unmapped complex: keep all interfaces (same rule as process_model)
            unmapped += 1
            writer.writerow(row)
            kept += 1
            continue
        chains = chain_set(row["native_chains"])
        if (chains & ab) and (chains & ag):
            # rewrite mapping columns to reference namespace (matches native_chains),
            # keeping the full per-complex sets as in the original format
            row["ab_chains"] = ",".join(sorted(ab))
            row["ag_chains"] = ",".join(sorted(ag))
            writer.writerow(row)
            kept += 1

print(f"{detailed.name}: {kept}/{total} interfaces kept "
      f"({unmapped} from unmapped complexes) -> {dst.name}")
