# Models Reference

## Protein Language Models (embedding generation)

### ESM-C 600M (`esmc_600m`)
- EvolutionaryScale, 36 layers, 1,152-dim embeddings
- Loaded via the `esm` SDK: `ESMC.from_pretrained("esmc_600m")`
- Embedding = mean-pool over residue positions of the final layer
- Used both directly AND as the input to ProteomeLM-S

### ProtT5-XL-UniRef50 (`Rostlab/prot_t5_xl_uniref50`)
- ProtTrans family, T5 encoder, 1,024-dim
- Requires space-separated amino acids; mean-pool excludes the EOS token
- Loaded with `T5EncoderModel` + `T5Tokenizer`

### ProteomeLM-S (`Bitbol-Lab/ProteomeLM-S`)
- 36.9M params, DistilBERT-style, 512-dim output
- **Takes ESM-C embeddings as input** (not raw sequence) — it contextualizes
  proteins against each other via attention
- The skill runs a **manual PyTorch forward pass** using the downloaded weights,
  because the `proteomelm` package's forward call is incompatible with
  transformers v5 (it passes `x=` where v5 expects `hidden_states=`, and calls a
  removed `get_head_mask`). The manual forward reproduces the architecture:
  `embedding_main + embedding_encoder` projection → 6 transformer blocks
  (pre-norm self-attention + GELU FFN).
- All proteins are passed together as one "proteome" so attention can capture
  inter-protein structure.

### ProteinBERT (`Rostlab/prot_bert`)
- BERT-based, 1,024-dim
- Loaded with `BertModel` + `BertTokenizer` (the slow tokenizer avoids a
  sentencepiece conversion error with the fast tokenizer)
- Requires space-separated amino acids; mean-pool excludes CLS/SEP

### Combined embedding
Concatenated in a fixed, deterministic order:
`ESM-C (1152) | ProtT5 (1024) | ProteomeLM (512) | ProteinBERT (1024) = 3,712`

## Classifiers

Nine models are trained and compared:

### Base (6)
| Model | Key settings |
|-------|--------------|
| Random Forest | 500 trees, max_depth 10, class_weight |
| Gradient Boosting | 300 est, depth 4, lr 0.05, subsample 0.8, sample_weight |
| SVM (RBF) | C=10, gamma scale, class_weight, probability |
| Logistic Regression | L2, lbfgs, class_weight |
| Extra Trees | 500 trees, max_depth 10, class_weight |
| MLP | (256,128), early stopping, sample_weight |

### Ensembles (3)
| Model | Description |
|-------|-------------|
| Soft Voting | Averages predicted probabilities of RF/SVM/LR/ExtraTrees/MLP |
| **Stacking (LR meta)** | LogReg meta-learner over the 5 base models via 5-fold internal CV. **Usually the best.** |
| Balanced Bagging | 20 ExtraTrees bags on balanced subsets (imbalanced-learn) |

## Why Stacking wins

The stacking meta-learner learns *how much to trust each base model* per region
of feature space, and the LogReg meta-learner produces well-calibrated
probabilities. This yields a moderate, stable optimal threshold and the best
precision/recall balance. Voting is a strong, simpler alternative. Tree-only and
bagging models tend to over-predict positives (high recall, low precision).

## Threshold per model

Each model's probability distribution differs, so the F1-optimal threshold is
tuned per model on cross-validation:
- SVM probabilities compress near 0/1 → low threshold (~0.1)
- MLP is overconfident → high threshold (~0.75)
- Stacking is well-calibrated → moderate threshold (~0.28)

The saved `trained_pipeline.pkl` stores the best model's threshold; predictions
apply it automatically.
