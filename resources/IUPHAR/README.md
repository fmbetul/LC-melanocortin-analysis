# IUPHAR Gene Lists

## Source

These files were downloaded from the IUPHAR/BPS Guide to Pharmacology (GtoPdb):  
https://www.guidetopharmacology.org/

**Download date:** October 2025  
**Database version:** GtoPdb 2024.3 (release December 2024)

## Files

| File | Contents | Used in |
|---|---|---|
| `GPCRTargets.csv` | All human GPCR gene symbols from the IUPHAR GPCR database | `02_meta_integration.Rmd`, Supp Figs 7 & 8 |
| `GtP_to_HGNC_mapping.csv` | IUPHAR target → HGNC symbol mapping (all receptor/ligand families) | `02_meta_integration.Rmd`, `04_visium_integrated.Rmd` |

## How to Re-download

1. Navigate to https://www.guidetopharmacology.org/download.jsp
2. Download "Targets and families" → export as CSV
3. For GPCR-specific list: filter by "GPCR" receptor family

Note: column names in newer releases may differ from those used here
(`HGNC.symbol` and `HGNC.Symbol`). Update the `read.csv` calls in the
analysis scripts if needed.