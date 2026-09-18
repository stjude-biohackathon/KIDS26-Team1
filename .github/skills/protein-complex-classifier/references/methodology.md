# Methodology & Design Rationale

This document captures the design decisions learned from the INO80/SRCAP
chromatin remodeling complex case study that motivated this skill.

## Dataset construction

### Positive set
Known complex interactors — curated subunits of the target complex(es). For
INO80/SRCAP: the ATPase core, actin-related proteins (ACTR5/8/6, ACTL6A),
RUVBL1/2, YEATS4, VPS72, ZNHIT1, etc. Start with ~20-50 high-confidence members.

### Negative set
Proteins known *not* to interact with the complex. Best practice: use nuclear
proteins from *other* complexes (spliceosome, nuclear pore, replication, DNA
repair, transcription machinery). This makes the classification task meaningful
— the model must distinguish complex membership, not just "nuclear vs
cytoplasmic". A cytoplasmic/metabolic negative set makes the task trivially easy
and the model won't generalize.

### Quality control
Before training, audit the gene sets:
- **Localization check**: confirm negatives are nuclear (bioservices/UniProt
  `cc_subcellular_location`). A negative in the wrong compartment weakens the
  contrast.
- **Gene symbol ambiguity**: some symbols map to unexpected UniProt entries
  (e.g. "ATR" can match a natriuretic peptide receptor instead of the kinase;
  "ALB" can match ALMS1 instead of albumin). Flag any gene whose fetched entry
  lacks expected annotations.
- **Overlap check**: ensure no gene appears in both positive and negative sets.

## Class imbalance (the central challenge)

Protein-complex datasets are inherently imbalanced — few known interactors vs
many possible non-interactors. In the case study, the ratio went from 1:10.6
(23 pos : 244 neg) to 1:4.6 (53 : 244) after adding more positives.

### Four simultaneous strategies
1. **Class weights** = `{0: 1.0, 1: neg/pos}` — the exact imbalance ratio,
   rounded to 1 decimal (e.g. 244/53 = 4.60 → 4.6). This was chosen over
   sklearn's `'balanced'` because it's interpretable and directly ties to the
   data. **Note:** the class weight interacts sensitively with threshold
   optimization — an unrounded 4.6038 vs rounded 4.6 can shift the CV-optimal
   threshold enough to change which model wins when several cluster near the top.
   The skill rounds to 1 decimal by default for reproducibility; override with
   `--class-weight <value>`.
2. **SMOTE** — synthetic minority oversampling, applied *inside each CV fold* and
   on the final training set, **never on test/validation**. Applied after PCA
   (works better in the reduced space).
3. **Stratified splits** — preserve the class ratio in train/test.
4. **Threshold optimization** — tune the decision boundary to maximize F1 rather
   than using a fixed 0.5. For imbalanced data the optimal threshold is often far
   from 0.5 (the case-study Stacking model used 0.28).

### What did NOT help
- Class weights alone (no SMOTE, no threshold tuning): recall stayed ~35%.
- A fixed 0.5 threshold with a tiny test set (5 positives): every model scored
  F1=0 on test because the threshold missed all positives even though ROC-AUC
  was 0.97. **Lesson: with few positives, report CV metrics and use threshold
  optimization; a single held-out split is unstable.**

## Dimensionality reduction

3,712 PLM features (or 3,860 with biological) vs ~200-300 samples → severe
overfitting risk. **PCA retaining 95% variance** reduces to ~120 components,
which:
- Eliminated overfitting (tree models were memorizing before PCA)
- Made SMOTE more effective (interpolation in a denser space)
- Sped up training

## Biological features (the key improvement)

PLM embeddings capture sequence-level signal but miss explicit functional
annotations. Adding 148 biological features improved the case-study Stacking
classifier:

| Metric | PLM only (3,712) | PLM + Bio (3,860) |
|--------|------------------|-------------------|
| Precision | 92.3% | **100%** |
| Recall | 75.0% | **81.2%** |
| F1 | 0.828 | **0.897** |
| MCC | 0.801 | **0.884** |

The most informative biological categories were **domains** (InterPro/Pfam
chromatin-domain flags), **GO terms** (chromatin remodeling, ATPase, histone
binding), and **PPI network** (shared STRING interactors with known positives).

## Model selection

Stacking (Logistic Regression meta-learner over RF/SVM/LogReg/ExtraTrees/MLP)
consistently won. It combines the diverse strengths of the base models and the
LR meta-learner calibrates the final probabilities well, giving a moderate,
stable optimal threshold. Soft Voting was a close second. Tree-only models
(RF/ExtraTrees) and Balanced Bagging over-predicted positives (low precision).

## When hyperparameter tuning won't help

If the misclassified positives have probabilities *far* below the threshold
(e.g. 0.02-0.06 when the boundary is 0.28, with a 0.85 gap to the nearest true
positive), those are genuinely hard cases — small accessory subunits whose
sequence/annotation profile resembles the negatives. Hyperparameter tuning shifts
probabilities marginally and won't close such a gap. Better levers: more positive
training examples, or co-complex/structural features.

## Reproducibility

- `random_state=42` throughout.
- All API results cached (`output/cache/`) — reruns are deterministic and fast.
- The trained pipeline bundles scaler + PCA + model + threshold so prediction on
  new genes reuses the exact training transforms.
