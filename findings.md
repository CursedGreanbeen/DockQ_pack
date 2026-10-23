# DockQ Package Findings

## Overview
The DockQ package is a batch quality assessment tool for AlphaFold predictions, specifically designed for antibody-antigen interface analysis.

## File Structure
```
/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/scripts/dockq/
├── cli.py          # Main CLI entry point
├── mapping.py      # AB/AG chain mapping utilities
├── analyze.py      # DockQ output parsing and analysis
├── run.py          # DockQ execution wrapper
├── resolution.py   # (not yet analyzed)
├── __init__.py     # Package init
└── __pycache__/    # Compiled Python
```

## Key Functions

### cli.py
- **`main()`**: Entry point that orchestrates the batch DockQ assessment
- **`find_model_files()`**: Discovers `*_model.cif` files in s3data directory
- **`find_reference()`**: Locates reference CIF files
- **`write_interface_report()`**: Writes AB-AG interfaces only
- **`append_results()`**: Appends results to TSV output

### mapping.py
- **`load_anarci_report()`**: Loads ANARCI crop report TSV to get chain type mappings
  - Status `CROPPED` → Antibody chain (AB)
  - Status `NOT_RECOGNIZED` → Antigen chain (AG)
  - Status `NO_ANTIBODY_CHAINS` → Skipped
- **`get_chain_type_mapping()`**: Returns `(ab_chains, ag_chains)` lists for a PDB ID
- **`build_ab_ag_pairs()`**: Creates all AB-AG chain pairs for DockQ --mapping flags
- **`format_dockq_mappings()`**: Formats pairs as `--mapping MODEL_CHAINS:NATIVE_CHAINS`

### analyze.py
- **`parse_dockq_output()`**: Parses DockQ stdout to extract:
  - Global DockQ score
  - Per-interface metrics: `native_chains`, `dockq_score`, `iRMSD`, `LRMSD`, `fnat`
- **`filter_ab_ag_interfaces()`**: Filters interfaces to keep only AB-AG interfaces
- **`is_ab_ag_interface()`**: Checks if interface contains both AB and AG chains
- **`determine_quality_bin()`**: Classifies DockQ score:
  - High: ≥ 0.80
  - Medium: ≥ 0.49
  - Acceptable: ≥ 0.23
  - Incorrect: < 0.23
- **`aggregate_interface_metrics()`**: Aggregates metrics by group (pdb_id, ab_chains, ab_ag_pair)

## DockQ Formula (Verified)
```
DockQ = (fnat + iRMSD_term + LRMSD_term) / 3

Where:
  iRMSD_term = 1 / (1 + (iRMSD/1.5)²)
  LRMSD_term = 1 / (1 + (LRMSD/8.5)²)
```

## Output Format
TSV with columns:
```
pdb_id | model_name | reference_path | interface_id | native_chains | ab_chains | ag_chains | dockq_score | iRMSD | LRMSD | fnat | quality_bin
```

## Critical Issue Found: Chain Mapping Problem

### Root Cause
The `ab_chains` and `ag_chains` columns show "N/A" because the ANARCI report path is not being found or the chain mapping is failing.

**Default ANARCI report path in cli.py:**
```python
ANARCI_REPORT_PATH = Path("/home/mullagaliamova/ClaudeWorkspace/PROJECTS/cdr-h3-folding/results/reports/anarci_crop_report.tsv")
```

**Problem:** The script is being run with a different ANARCI report path (or the file doesn't exist), so `ab_ag_mapping` is empty, resulting in `ab_chains` and `ag_chains` being "N/A".

### Interface Filtering Logic
The code correctly identifies AB-AG interfaces using `is_ab_ag_interface()`:
- An interface must contain **at least one AB chain AND at least one AG chain**
- Interfaces with only AB chains (e.g., heavy-light chain interface B-C) are NOT AB-AG interfaces

### For 8iuz Specifically
**Chain composition (from FASTA):**
- Chain A: Antigen (240 residues)
- Chain B: Antibody heavy chain (118 residues) - VH domain
- Chain C: Antibody light chain (107 residues) - VL domain

**Expected mappings:**
- AB chains: B, C
- AG chains: A

**Interfaces reported:**
| Interface | native_chains | Expected Type |
|-----------|---------------|---------------|
| iface1 | A, B | AB-AG ✓ |
| iface2 | A, C | AB-AG ✓ |
| iface3 | B, C | Intra-AB (heavy-light) ✗ |

**Problem:** iface3 (B, C) is the interface between antibody heavy and light chains. This is NOT an antibody-antigen interface and should be:
1. Excluded from the AB-AG interface report, OR
2. Flagged as intra-antibody

The filtering logic in `cli.py` lines 428-430 does filter for `is_ab_ag`:
```python
if args.interface_report:
    ab_ag_results = [r for r in results if r.get("is_ab_ag", False)]
    append_results(ab_ag_results, args.interface_report)
```

But this only works if `ab_chains` and `ag_chains` are properly populated from the ANARCI report.

## Recommendations

1. **Verify ANARCI report path** - Ensure the `--anarci-report` flag points to the correct file
2. **Check ANARCI report content** - Verify 8iuz entries exist with correct chain mappings
3. **Consider fallback logic** - If no ANARCI mapping, could use chain sequence analysis to identify antibody chains (VH/VL signature sequences)
4. **Add validation** - Warn if `ab_chains` is empty for a complex that should have antibodies

## DockQ Calculation Verification
Tested 5 random interfaces - all DockQ scores calculated correctly using the standard formula.
