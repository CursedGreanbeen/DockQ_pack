#!/usr/bin/env python3
"""Recreate an AB-AG interface-only report from a full per-interface detailed report.

Usage: recreate_interface_report.py <detailed_report.tsv> <interface_report.tsv>
"""
import csv, sys
from pathlib import Path

sys.path.insert(0, "/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/scripts")  # parent of the dockq package
from dockq.analyze import is_ab_ag_interface

src = Path(sys.argv[1])   # full per-interface report
dst = Path(sys.argv[2])   # recreated AB-AG report

def split(s):
    return [c for c in s.split(",") if c and c != "N/A"]

kept = total = 0
with open(src, newline="") as fi, open(dst, "w", newline="") as fo:
    reader = csv.DictReader(fi, delimiter="\t")
    writer = csv.DictWriter(fo, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in reader:
        total += 1
        ab, ag = split(row["ab_chains"]), split(row["ag_chains"])
        # same logic as process_model: unmapped complexes keep all interfaces
        if not (ab and ag) or is_ab_ag_interface(row["native_chains"], ab, ag):
            writer.writerow(row)
            kept += 1

print(f"{src.name}: {kept}/{total} interfaces written to {dst.name}")
