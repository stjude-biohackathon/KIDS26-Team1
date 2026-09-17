# Raw DepMap inputs (Steps 1 & 2 source)

These are the unmodified DepMap files the pipeline reads. The two large matrices
feed `analysis/build_features.py`; `Model.csv` is the sample-metadata table that
supplies every cell line's cancer type and patient age.

| file | rows × cols | content |
|---|---|---|
| `CRISPRGeneEffect.csv` | 1208 cell lines × 18,531 genes | Chronos gene-effect scores (more negative = more essential) |
| `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | 1775 profiles × 19,221 genes | log2(TPM+1) expression, protein-coding genes (profile-level; collapsed to one row per `ModelID` during harmonization) |
| `Model.csv` | 2154 cell lines × 49 fields | Sample metadata — `ModelID`, `OncotreeLineage`, `OncotreePrimaryDisease`, `AgeCategory` (Pediatric/Adult/Fetus/Unknown), etc. |

Column headers on the two matrices use DepMap's `SYMBOL (EntrezID)` convention;
`build_features.py` strips the Entrez suffix during harmonization.

`Model.csv` is copied verbatim (no filtering/renaming/derivation) into
`data/depmap_processed/model_meta.parquet` — the only change is CSV → Parquet
serialization; all 2,154 × 49 cells are identical. That parquet drives the
pediatric/adult split and the cohort cancer-type figures
(`analysis/cohort_summary.py`).

**Source:** DepMap 26Q1 Public, downloaded manually from the Broad Institute
DepMap portal (https://depmap.org/portal/data_page/?tab=allData). File
dimensions and first-50 MB MD5 checksums are recorded in
`data/analysis_outputs/data_manifest.csv`. Please cite DepMap when using these data.

**Regenerate the feature matrix:**
```bash
cd <package root>
python analysis/build_features.py \
    --gene-effect data/depmap_raw/CRISPRGeneEffect.csv \
    --expression  data/depmap_raw/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv \
    --seeds       data/training/seed_set.csv \
    --out         data/depmap_processed/features.parquet
```
This yields 16 features over 19,283 genes (`features.parquet`).
