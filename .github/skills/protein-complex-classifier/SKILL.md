---
name: protein-complex-classifier
description: End-to-end pipeline that predicts whether proteins interact with a target protein complex (e.g. INO80/SRCAP chromatin remodeling complexes) from gene lists. Takes positive (known interactors) and negative (non-interactors) gene symbol lists, fetches sequences from UniProt, generates protein language model embeddings (ESM-C 600M, ProtT5-XL, ProteomeLM-S, ProteinBERT), extracts biological features (physicochemical, sequence composition, InterPro domains, GO terms, STRING PPI network), builds a combined feature matrix, and trains an ensemble ML classifier (Stacking/Voting/Bagging with PCA + SMOTE + class weighting + threshold optimization). Use when the user wants to build a protein-complex interaction predictor, classify proteins as complex interactors/non-interactors, generate protein embeddings from gene lists, or run predictions on new gene lists with a trained model. Trigger on requests involving protein-complex interaction prediction, protein embeddings from gene symbols, or chromatin remodeling complex membership classification.
license: MIT
compatibility: Requires Python 3.11-3.13 with uv, ~8GB RAM (16GB recommended), internet access to UniProt/InterPro/QuickGO/STRING/HuggingFace. CPU-only OK (~30-60 min full run for ~300 proteins). First run downloads ~6GB of model weights. Reproducing the deployed model (see below) needs no internet access at all.
metadata:
  version: "1.0"
  author: Generated from INO80/SRCAP classifier workflow
---

# Protein-Complex Interaction Classifier

End-to-end workflow that predicts whether proteins interact with a target protein
complex, starting from gene symbol lists and ending with a trained ensemble classifier
and predictions on new genes.

## What this skill does

```
Gene lists (positive + negative)
    │
    ├─▶ 1. Fetch sequences from UniProt
    │
    ├─▶ 2. Generate PLM embeddings (4 models → 3,712 dims)
    │       ESM-C 600M · ProtT5-XL · ProteomeLM-S · ProteinBERT
    │
    ├─▶ 3. Extract biological features (148 dims)
    │       physicochemical · composition · domains · GO · PPI
    │
    ├─▶ 4. Build combined feature matrix (3,860 dims)
    │
    ├─▶ 5. Train ensemble classifier
    │       PCA(95%) + SMOTE + class_weight + threshold optimization
    │       9 models: RF, GBM, SVM, LogReg, ExtraTrees, MLP,
    │                 Soft Voting, Stacking, Balanced Bagging
    │
    └─▶ 6. Predict on new (unlabeled) gene lists
```

## Setup (once)

```bash
cd <skill_dir>
bash scripts/setup.sh
```

The skill ships its own `pyproject.toml` + `uv.lock` (pinned to the exact
versions validated for this pipeline: numpy 2.3.5, torch 2.11.0, transformers
4.57.6, esm 3.4.0, etc.). `setup.sh` runs `uv sync` against this lock file —
not an open-ended `uv add` — so a fresh install reproduces the same
environment rather than resolving to whatever's newest. The `esm` package is
installed separately with `--no-deps` (torchtext/torchvision lack Python 3.13
wheels and are not needed for ESM-C embeddings). See
[references/dependencies.md](references/dependencies.md).

## Input format

Two plain-text files, one gene symbol per line (HGNC symbols, human):

```
# positive_set                # negative_set
INO80                         SMARCA1
SRCAP                         BAZ1A
RUVBL1                        CHD1
...                           ...
```

Gene symbols are mapped to reviewed human UniProt entries automatically.

## Usage

### Reproducing the deployed INO80/SRCAP model exactly (frozen features, no live APIs)

**Use this deployed model**
`data/models/trained_pipeline.pkl` whenever you want to verify a result, 
to audit the classifier, or to build predictions from the exact production pipeline.
Do **not** use this path if you want to retrain on an updated/expanded gene
set or deliberately pick up newer GO/STRING annotations — use the "Full
training pipeline" section below for that.

**Why this path exists:** the biological features (GO terms, InterPro
domains, STRING PPI scores) are fetched live from public databases that
change over time. Re-running the pipeline from gene lists at a later date
recomputes ~32 of the 3,860 feature columns from whatever UniProt/STRING
say *today*, which is enough to shift the CV-optimal decision threshold and
a couple of test-set predictions even though the classifier code itself is
perfectly deterministic (verified: given an identical feature matrix, this
skill's `train_classifier.py` reproduces `trained_pipeline.pkl`'s test-set
predictions bit-for-bit — 0.0 max probability difference across all 90
held-out genes). The fix isn't "try harder to match live data" — it's to
never re-fetch it when exact reproduction is the goal. This skill ships a
snapshot of the exact data that produced the deployed model:

```
data/
├── positive_set                                 # 53 genes (dataset composition, for reference)
├── negative_set                                  # 244 genes (dataset composition, for reference)
├── embeddings/combined_feature_matrix.npz        # 297 x 3,712 PLM embeddings
├── features/full_feature_matrix.npz              # 297 x 3,860 full feature matrix (PLM + bio)
├── features/full_feature_matrix.csv              # same, human-readable (needed by predict.py --train-features)
├── features/_string_cache.json                   # training-time STRING partners cache (needed by predict.py --train-string-cache)
└── models/trained_pipeline.pkl             # the actual deployed classifier — use directly if you don't need to retrain
```

**Option A — use the deployed model directly (no training at all):**

```bash
uv run python scripts/predict.py \
    --genes <your_gene_list> \
    --model data/models/trained_pipeline.pkl \
    --train-features data/features/full_feature_matrix.csv \
    --train-embeddings data/embeddings/combined_feature_matrix.npz \
    --train-string-cache data/features/_string_cache.json \
    --output predictions.csv
```

This is the fastest and most direct way to get predictions from the exact
model the pickle represents; it still fetches sequences/embeddings/features
live for your *new* query genes (that's unavoidable — they aren't in the
frozen set) but the training side is fully pinned.

**Option B — retrain from the frozen features (to confirm reproducibility, inspect the training code, or regenerate `all_trained_models.pkl`):**

```bash
mkdir -p output/features output/embeddings
cp data/features/full_feature_matrix.npz output/features/
cp data/features/full_feature_matrix.csv output/features/
cp data/features/_string_cache.json output/features/
cp data/embeddings/combined_feature_matrix.npz output/embeddings/

uv run python scripts/train_classifier.py \
    --output-dir output --split 0.3 --pca-variance 0.95 \
    --class-weight auto --best-metric F1-Score
```

Expect: 297 genes (53 pos / 244 neg) → 3,860 features → PCA to 120 components →
best model **Stacking (LR meta)** at F1 = 0.897 (threshold = 0.28), Precision = 1.000,
Recall = 0.812 on the held-out test split — verified to reproduce exactly, run after
run, because no live data enters the computation.

### Full training pipeline (generic — train on a *different* complex, or deliberately refresh GO/STRING data)

Use this path to train on a new protein complex (no frozen data exists for
it), or if you specifically want the classifier to reflect the current state
of GO/InterPro/STRING rather than the frozen snapshot above. This hits live
UniProt/InterPro/STRING APIs and results will drift slightly from the
deployed model run to run (see Notes & caveats).

```bash
uv run python scripts/run_pipeline.py \
    --positive data/positive_set \
    --negative data/negative_set \
    --output-dir output/ \
    --complex-name "INO80/SRCAP"
```

To train on a genuinely different complex, replace the gene list paths and
`--complex-name` with your own:

```bash
uv run python scripts/run_pipeline.py \
    --positive <your_positive_gene_list> \
    --negative <your_negative_gene_list> \
    --output-dir output/ \
    --complex-name "<Your Complex Name>"
```

This produces:
- `output/embeddings/` — per-model + combined embedding CSV/NPZ
- `output/features/` — biological features + full feature matrix
- `output/models/trained_pipeline.pkl` — the production classifier
- `output/models/all_trained_models.pkl` — all 9 models for comparison
- `output/results/` — metrics CSVs + per-protein predictions
- `output/figures/` — confusion matrix, ROC, PR curves, comparison plots

### Run individual steps

```bash
# Step 1-2: sequences + embeddings only
uv run python scripts/generate_embeddings.py \
    --positive data/positive_set --negative data/negative_set \
    --output-dir output/

# Step 3-4: biological features (needs embeddings from step 2)
uv run python scripts/extract_biological_features.py --output-dir output/

# Step 5: train classifier (needs full_feature_matrix from step 4, or the
# frozen one copied in per "Option B" above)
uv run python scripts/train_classifier.py --output-dir output/ \
    --split 0.3 --pca-variance 0.95

# Step 6: predict on new genes with the trained model
uv run python scripts/predict.py \
    --genes data/validation_genes \
    --model output/models/trained_pipeline.pkl \
    --train-features output/features/full_feature_matrix.csv \
    --train-embeddings output/embeddings/combined_feature_matrix.npz \
    --output output/results/predictions.csv
```

`--train-features` and `--train-embeddings` are required for full (non-`--plm-only`)
models: they align new genes' 3,860 columns to the exact training order and identify
which training genes were positives (for PPI "linking to known positives" features).
`predict.py` auto-detects the matching STRING cache from `--train-embeddings`'
location; pass `--train-string-cache <path>` explicitly if it isn't found (a printed
warning tells you when this happens and which paths were checked) — this is exactly
why the frozen bundle above ships `features/_string_cache.json` alongside the feature
matrix.

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--split` | 0.3 | Test set fraction (stratified) |
| `--pca-variance` | 0.95 | Variance retained by PCA |
| `--class-weight` | auto | Positive class weight; `auto` = neg/pos ratio |
| `--use-smote` | true | SMOTE oversampling of minority class in training folds |
| `--optimize-threshold` | true | Tune decision threshold to maximize F1 |
| `--best-metric` | `F1-Score` | Metric for selecting the best model |

The class weight defaults to the exact imbalance ratio (e.g. 244 neg / 53 pos ≈ 4.6,
rounded to 1 decimal place — see Notes & caveats for why the rounding matters).
See [references/methodology.md](references/methodology.md) for the rationale behind
each design choice, learned from the INO80/SRCAP case study.

## Models

**Protein language models** (embeddings):

| Model | HF ID | Dim |
|-------|-------|-----|
| ESM-C 600M | `esmc_600m` (EvolutionaryScale) | 1,152 |
| ProtT5-XL | `Rostlab/prot_t5_xl_uniref50` | 1,024 |
| ProteomeLM-S | `Bitbol-Lab/ProteomeLM-S` | 512 |
| ProteinBERT | `Rostlab/prot_bert` | 1,024 |

**Classifiers** (ensemble): Random Forest, Gradient Boosting, SVM-RBF, Logistic
Regression, Extra Trees, MLP, Soft Voting, **Stacking (LR meta)**, Balanced Bagging.
Stacking typically wins. See [references/models.md](references/models.md).

## Biological features (148)

| Category | Count | Source |
|----------|-------|--------|
| Physicochemical | 33 | Sequence (BioPython ProtParam) |
| Sequence composition | 53 | Sequence (dipeptides, disorder, grouped AA) |
| Domain/motif | 25 | InterPro/Pfam API |
| GO terms | 28 | UniProt/QuickGO |
| PPI network | 9 | STRING database |

Adding biological features to PLM embeddings improved the INO80/SRCAP classifier from
F1=0.828 to F1=0.897 (precision 92%→100%, recall 75%→81%). See
[references/methodology.md](references/methodology.md).

## Handling class imbalance

Protein-complex datasets are typically imbalanced (few known interactors). This skill
addresses it four ways simultaneously:
1. **Class weights** — positive class weighted by the neg/pos ratio
2. **SMOTE** — synthetic minority oversampling (training folds only)
3. **Stratified splits** — preserve class ratio in train/test
4. **Threshold optimization** — tune the decision boundary, not fixed 0.5

Never apply SMOTE or resampling to test/validation data — the scripts enforce this.

## Notes & caveats

- All API results (sequences, InterPro, GO, STRING) are cached; reruns are fast.
- **GO/InterPro/STRING data drifts over time.** These are fetched live from public
  databases, so re-running the pipeline months later may yield slightly different
  values for the ~62 GO/domain/PPI features if UniProt/STRING have updated their
  annotations — this is expected and not a code bug. Physicochemical and
  sequence-composition features are deterministic (computed purely from the fetched
  sequence) and should reproduce exactly.
- Gene symbols must be current HGNC symbols; ambiguous symbols may map to the wrong
  UniProt entry (validate flagged genes — see [references/methodology.md](references/methodology.md)).
- The trained pipeline bundles scaler + PCA + model + threshold, so prediction on new
  genes reuses the exact training transforms.
- `predict.py`'s PPI "linking to known positives" features need the training-time
  STRING partner data for the *positive* genes (not just the query genes). Pass
  `--train-string-cache` explicitly if auto-detection doesn't find it (the script
  prints a warning listing the paths it checked); otherwise those 2 of 148 features
  silently default to 0.
- If you run `uv sync` after `setup.sh`, it will remove the out-of-band-installed
  `esm` package (it isn't a declared `pyproject.toml` dependency — see
  [references/dependencies.md](references/dependencies.md) for why). Re-run
  `uv pip install esm --no-deps` before generating embeddings again.
- For localization/QC auditing of gene sets, see
  [references/methodology.md](references/methodology.md#quality-control).
