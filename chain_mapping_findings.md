# DockQ Chain Mapping Findings

## Summary

When running DockQ on AlphaFold3 model predictions vs reference structures, **chain-based filtering of AB-AG interfaces is not working** because of chain annotation mismatches.

---

## Problem: No Filtering Between Reports

The detailed report and interface report are **identical** (same line count, same content) because:

1. **No AG chains in ANARCI report**: The reference ANARCI report (`results/reports/anarci_report.tsv`) contains:
   - 212 AB (antibody) chains
   - **0 AG (antigen) chains**
   
2. **Status codes not recognized**: The ANARCI report uses status `NO_ANTIBODY_DOMAIN` for antigen chains, but the mapping code only recognized:
   - `ANTIBODY` → AB
   - `NOT_RECOGNIZED` → AG (not present in our report!)

3. **Chain ID mismatch**: Even if we fix the status code issue:
   - Reference ANARCI chains for 11hk: A, D, E, F (E, F = antibody)
   - DockQ interface chains for 11hk: G, E, F
   - Chain G exists in models but not in reference ANARCI report

---

## Root Cause: Models Need Their Own ANARCI Annotation

The AlphaFold3 **models** have different chain naming than the **reference structures**:
- Models use chains like: A, B, C, D, E, F, G...
- References use different chain IDs

To properly filter AB-AG interfaces, we need to:
1. Extract FASTA sequences from each model CIF file
2. Run ANARCI on the model sequences
3. Generate a model-specific ANARCI report

---

## Technical Challenge: CIF to FASTA Extraction

AlphaFold3 outputs are in CIF format (not PDB). Extracting sequences requires:

### Option 1: Parse _atom_site records manually
- Read residue names from atom records
- Convert 3-letter codes to 1-letter
- Group by chain ID

### Option 2: Use gemmi library
```python
import gemmi
structure = gemmi.read_structure("model.cif")
for model in structure:
    for chain in model:
        for residue in chain:
            # residue.name gives 3-letter code (e.g., "ALA")
            # Need manual conversion to 1-letter
```

**Issue**: gemmi's `one_letter_code()` function expects a sequence list, not a single residue name. Simple dictionary lookup works better.

---

## Required Next Steps

1. **Create FASTA extraction script**: `/mnt/1857e392-c689-44a8-9c3e-e7f9beab82f6/scripts/extract_model_fasta.py`
   - Parse CIF files using gemmi
   - Extract per-chain sequences
   - Output FASTA format

2. **Run ANARCI on model FASTAs**:
   - Process all 10 PDB complexes in single-seed folder
   - Generate `anarci_model_report.tsv` with same format as reference report

3. **Update DockQ CLI** (optional):
   - Already fixed to recognize `NO_ANTIBODY_DOMAIN` as antigen
   - May need to handle chain ID mapping between models and references

---

## Files Created/Modified

| File | Purpose |
|------|---------|
| `scripts/extract_model_fasta.py` | Extract FASTA from AF3 CIF models |
| `scripts/dockq/mapping.py` | Fixed to recognize NO_ANTIBODY_DOMAIN as AG |
| `scripts/dockq/chain_mapping_findings.md` | This document |

---

## Current Status

- ✅ DockQ batch run completed (50 models, 160 interfaces)
- ✅ Quality assessment done (23% High, 22% Medium, 5% Acceptable, 50% Incorrect)
- ❌ Chain-based filtering NOT working (reports identical)
- ⏳ Need model ANARCI annotation to enable filtering
