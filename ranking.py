#!/usr/bin/env python3
"""
Model ranking / pre-filtering for the DockQ batch pipeline.

Keeps only the top N model samples per seed, ranked by ranking_score from
the per-complex file <pdb>_ranking_scores.csv (columns: seed,sample,
ranking_score). E.g. --top-per-seed 1 cuts 500 seeds x 5 samples down to
500 DockQ runs per complex.

ranking_score is the AlphaFold composite (0.8*ipTM + 0.2*pTM, plus
disorder/clash penalties in AF3-style pipelines) -- not pLDDT.
"""

import csv
import re
from pathlib import Path
from typing import Optional

# seed-1002927834_sample-0
_SEED_DIR_RE = re.compile(r'^seed-(?P<seed>.+)_sample-(?P<sample>\d+)$')


def parse_seed_sample(dir_name: str) -> Optional[tuple[str, int]]:
    """Parse (seed, sample) from a seed dir name like 'seed-1002927834_sample-0'."""
    m = _SEED_DIR_RE.match(dir_name)
    if not m:
        return None
    return m.group('seed'), int(m.group('sample'))


def load_top_samples_per_seed(ranking_file: Path, top_per_seed: int) -> Optional[dict[str, set[int]]]:
    """
    Read the ranking CSV once and return {seed: set of top-N sample ids}.

    Returns None if ranking_file does not exist.
    """
    if not ranking_file.exists():
        return None

    by_seed: dict[str, list[tuple[int, float]]] = {}
    with open(ranking_file) as f:
        for row in csv.DictReader(f):
            seed, sample, score = row.get('seed'), row.get('sample'), row.get('ranking_score')
            if not seed or sample is None or score is None:
                continue
            try:
                by_seed.setdefault(seed, []).append((int(sample), float(score)))
            except ValueError:
                continue

    return {
        seed: {s for s, _ in sorted(entries, key=lambda x: -x[1])[:top_per_seed]}
        for seed, entries in by_seed.items()
    }


def filter_models_by_ranking(models: list[dict], top_per_seed: int) -> list[dict]:
    """
    Keep only the top `top_per_seed` models per seed (by ranking_score).

    Each PDB's ranking CSV is read once and cached. Models are kept with a
    warning -- never silently dropped -- when their seed/sample can't be
    parsed, their seed is missing from the CSV, or the ranking file is
    absent.
    """
    cache: dict[Path, Optional[dict[str, set[int]]]] = {}
    kept = []

    for model in models:
        model_path = Path(model['model_path'])
        parsed = parse_seed_sample(model_path.parent.name)

        if parsed is None:
            print(f"WARNING: cannot parse seed/sample from {model_path.name}; keeping")
            kept.append(model)
            continue

        seed, sample = parsed
        # <task>/<pdb>/seed-<seed>_sample-<k>/<file> -> ranking CSV sits in <pdb>/
        ranking_file = model_path.parent.parent / f"{model['pdb_id']}_ranking_scores.csv"

        if ranking_file not in cache:
            cache[ranking_file] = load_top_samples_per_seed(ranking_file, top_per_seed)
            if cache[ranking_file] is None:
                print(f"WARNING: ranking file not found: {ranking_file}; no filtering for this PDB")

        top = cache[ranking_file]
        if top is None:
            kept.append(model)
        elif seed not in top:
            print(f"WARNING: {model['pdb_id']}: seed {seed} missing from ranking file; keeping")
            kept.append(model)
        elif sample in top[seed]:
            kept.append(model)

    return kept
