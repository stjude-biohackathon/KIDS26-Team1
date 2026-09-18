# Dependencies

## Installation

Run `bash scripts/setup.sh` once. It installs everything into a uv-managed
environment.

## Why `esm` needs `--no-deps`

The `esm` package (EvolutionaryScale, for ESM-C) declares `torchtext` and
`torchvision` as hard dependencies. These lack wheels for Python 3.13 and are
**not used** by ESM-C embedding generation. Installing `esm --no-deps` and
providing its actual runtime deps explicitly (einops, biotite, brotli, zstd,
tenacity, cloudpathlib, msgpack-numpy, httpx) avoids the build failure.

## Full dependency list

| Package | Purpose |
|---------|---------|
| numpy, pandas, scipy | Core data |
| scikit-learn | Classifiers, PCA, metrics |
| imbalanced-learn | SMOTE, BalancedBagging |
| matplotlib, seaborn | Figures |
| torch | Model inference backend |
| transformers | ProtT5, ProteinBERT, tokenizers |
| huggingface-hub, safetensors | ProteomeLM weight download |
| accelerate | ESM-C model loading |
| sentencepiece, tokenizers, protobuf | ProtT5/BERT tokenizers |
| esm (`--no-deps`) | ESM-C 600M |
| einops, biotite, brotli, zstd, tenacity, cloudpathlib, msgpack-numpy, httpx | esm runtime deps |
| biopython | ProtParam physicochemical features |
| bioservices | UniProt/GO queries |
| requests, joblib | HTTP + model serialization |

## System requirements

- Python 3.11-3.13
- ~8 GB RAM (16 GB recommended for ProtT5-XL)
- ~6 GB disk for cached model weights (`~/.cache/huggingface/`)
- Internet access to: rest.uniprot.org, ebi.ac.uk (InterPro), string-db.org,
  huggingface.co
- CPU-only is fine (full ~300-protein run ≈ 30-60 min; first run adds model
  download time)

## Model weight sources (downloaded on first use)

| Model | Source | Approx size |
|-------|--------|-------------|
| ESM-C 600M | HuggingFace `esmc_600m` | ~2.3 GB |
| ProtT5-XL | `Rostlab/prot_t5_xl_uniref50` | ~2.4 GB |
| ProteomeLM-S | `Bitbol-Lab/ProteomeLM-S` | ~74 MB |
| ProteinBERT | `Rostlab/prot_bert` | ~1.6 GB |
