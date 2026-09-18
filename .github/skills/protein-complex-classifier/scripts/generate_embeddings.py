#!/usr/bin/env python3
"""
Step 1-2: Fetch sequences + generate PLM embeddings.
Produces per-model and combined embedding CSV/NPZ files.
"""
import os, sys, json, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline_lib as lib


def main():
    ap = argparse.ArgumentParser(description="Generate PLM embeddings from gene lists")
    ap.add_argument("--positive", required=True, help="Positive gene list file")
    ap.add_argument("--negative", required=True, help="Negative gene list file")
    ap.add_argument("--output-dir", default="output", help="Output directory")
    args = ap.parse_args()

    embed_dir = os.path.join(args.output_dir, "embeddings")
    cache_dir = os.path.join(args.output_dir, "cache")
    os.makedirs(embed_dir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)

    # ── Load genes ──
    pos_genes = lib.read_genes(args.positive)
    neg_genes = lib.read_genes(args.negative)
    all_genes = pos_genes + neg_genes
    labels = [1] * len(pos_genes) + [0] * len(neg_genes)
    print(f"Genes: {len(pos_genes)} positive + {len(neg_genes)} negative = {len(all_genes)}")

    # ── Fetch sequences ──
    print("\nFetching sequences from UniProt...")
    seq_data = lib.fetch_sequences(all_genes,
                                   cache_path=os.path.join(cache_dir, "sequences.json"))

    gene_names, seqs, y, uid_map = [], [], [], {}
    for gene, label in zip(all_genes, labels):
        if gene in seq_data and seq_data[gene].get("sequence"):
            gene_names.append(gene); seqs.append(seq_data[gene]["sequence"])
            y.append(label); uid_map[gene] = seq_data[gene]["uniprot_id"]
    y = np.array(y)
    N = len(gene_names)
    print(f"Ready: {N} proteins ({int(y.sum())} pos, {N - int(y.sum())} neg)")

    def save(name, emb, dim):
        cols = [f"dim_{j}" for j in range(dim)]
        df = pd.DataFrame(emb, columns=cols)
        df.insert(0, "Gene", gene_names)
        df.insert(1, "Label", ["positive" if l == 1 else "negative" for l in y])
        df.insert(2, "UniProt_ID", [uid_map[g] for g in gene_names])
        safe = name.lower().replace(" ", "_").replace("-", "_")
        df.to_csv(os.path.join(embed_dir, f"{safe}_embeddings.csv"), index=False)
        np.savez(os.path.join(embed_dir, f"{safe}_embeddings.npz"),
                 X=emb, y=y, genes=gene_names)
        print(f"  -> {safe}_embeddings.csv ({N} x {dim})")

    # ── 4 models ──
    print("\n[1/4] ESM-C 600M...")
    esmc, esmc_dim = lib.embed_esmc(seqs, gene_names)
    save("esmc_600m", esmc, esmc_dim)

    print("\n[2/4] ProtT5-XL...")
    t5, t5_dim = lib.embed_prott5(seqs, gene_names)
    save("prottrans_prott5_xl", t5, t5_dim)

    print("\n[3/4] ProteomeLM-S...")
    plm, plm_dim = lib.embed_proteomelm(esmc)
    save("proteomelm_s", plm, plm_dim)

    print("\n[4/4] ProteinBERT...")
    bert, bert_dim = lib.embed_proteinbert(seqs, gene_names)
    save("proteinbert", bert, bert_dim)

    # ── Combined ──
    print("\nBuilding combined feature matrix...")
    combined = np.concatenate([esmc, t5, plm, bert], axis=1)
    order = [("esmc_600m", esmc_dim), ("prottrans_prott5_xl", t5_dim),
             ("proteomelm_s", plm_dim), ("proteinbert", bert_dim)]
    cols = []
    for name, d in order:
        cols += [f"{name}_dim_{j}" for j in range(d)]
    cdf = pd.DataFrame(combined, columns=cols)
    cdf.insert(0, "Gene", gene_names)
    cdf.insert(1, "Label", ["positive" if l == 1 else "negative" for l in y])
    cdf.insert(2, "UniProt_ID", [uid_map[g] for g in gene_names])
    cdf.to_csv(os.path.join(embed_dir, "combined_feature_matrix.csv"), index=False)
    np.savez(os.path.join(embed_dir, "combined_feature_matrix.npz"),
             X=combined, y=y, genes=gene_names,
             model_order=[n for n, _ in order], dims=[d for _, d in order])
    print(f"  -> combined_feature_matrix.csv ({N} x {combined.shape[1]})")
    print(f"\nDone. Embeddings in {embed_dir}/")


if __name__ == "__main__":
    main()
