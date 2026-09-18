#!/usr/bin/env python3
"""
Step 6: Predict on a new (unlabeled) gene list with a trained pipeline.
Generates embeddings + biological features, then applies the saved model.
"""
import os, sys, json, time, argparse
import numpy as np
import pandas as pd
import joblib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline_lib as lib


def main():
    ap = argparse.ArgumentParser(description="Predict interactions for new genes")
    ap.add_argument("--genes", required=True, help="Gene list file (unlabeled)")
    ap.add_argument("--model", required=True, help="trained_pipeline.pkl")
    ap.add_argument("--output", default="predictions.csv")
    ap.add_argument("--train-features", help="Training full_feature_matrix.csv (for column alignment). "
                    "Required unless --plm-only model.")
    ap.add_argument("--train-embeddings", help="Training combined_feature_matrix.npz (identifies which "
                    "training genes were positives, for the PPI 'shared/linking-to-positives' features).")
    ap.add_argument("--train-string-cache", default=None,
                    help="Path to the training-time STRING-partners cache JSON (a dict of "
                    "gene -> list of {partner, score}) for the KNOWN POSITIVE genes. Needed to "
                    "compute 'n_positive_genes_linking_to_this' / 'frac_positives_linking' correctly. "
                    "If omitted, common locations are auto-detected relative to --train-embeddings "
                    "(<output-dir>/cache/string.json from this skill's own pipeline, or "
                    "<data-dir>/features/_string_cache.json from the original project scripts). "
                    "If none is found, those 2 of 148 features default to 0 and a warning is printed.")
    ap.add_argument("--cache-dir", default="predict_cache")
    args = ap.parse_args()

    os.makedirs(args.cache_dir, exist_ok=True)
    pipeline = joblib.load(args.model)
    plm_only = pipeline.get("plm_only", False)
    print(f"Model: {pipeline['best_model_name']} (threshold={pipeline['optimal_threshold']}, "
          f"plm_only={plm_only})")

    # ── Genes + sequences ──
    genes = lib.read_genes(args.genes)
    print(f"\nFetching sequences for {len(genes)} genes...")
    seq_data = lib.fetch_sequences(genes, cache_path=os.path.join(args.cache_dir, "sequences.json"))
    gene_names = [g for g in genes if g in seq_data and seq_data[g].get("sequence")]
    seqs = [seq_data[g]["sequence"] for g in gene_names]
    uid_map = {g: seq_data[g]["uniprot_id"] for g in gene_names}
    N = len(gene_names)
    print(f"Ready: {N} proteins")

    # ── PLM embeddings ──
    print("\nGenerating PLM embeddings...")
    print("  [1/4] ESM-C..."); esmc, ed = lib.embed_esmc(seqs, gene_names, verbose=False)
    print("  [2/4] ProtT5..."); t5, td = lib.embed_prott5(seqs, gene_names, verbose=False)
    print("  [3/4] ProteomeLM..."); plm, pd_ = lib.embed_proteomelm(esmc)
    print("  [4/4] ProteinBERT..."); bert, bd = lib.embed_proteinbert(seqs, gene_names, verbose=False)
    X_plm = np.concatenate([esmc, t5, plm, bert], axis=1)

    plm_cols = []
    for name, d in [("esmc_600m", ed), ("prottrans_prott5_xl", td),
                    ("proteomelm_s", pd_), ("proteinbert", bd)]:
        plm_cols += [f"{name}_dim_{j}" for j in range(d)]
    plm_df = pd.DataFrame(X_plm, columns=plm_cols)

    # ── Biological features (unless plm_only) ──
    if plm_only:
        val = plm_df
    else:
        print("\nExtracting biological features...")
        from bioservices import UniProt as UniProtService
        u = UniProtService(verbose=False)

        # Known positives + their STRING partners (for PPI features)
        known_pos, train_string = set(), {}
        if args.train_embeddings and os.path.exists(args.train_embeddings):
            tnpz = np.load(args.train_embeddings, allow_pickle=True)
            tg, ty = list(tnpz["genes"]), tnpz["y"]
            known_pos = set(tg[i] for i in range(len(tg)) if ty[i] == 1)

        # Locate the training-time STRING cache for the known positives.
        # Explicit --train-string-cache always wins; otherwise auto-detect
        # the two layouts this pipeline is known to produce/consume:
        #   1. <output-dir>/cache/string.json      (this skill's own run_pipeline.py output)
        #   2. <data-dir>/features/_string_cache.json  (the original project scripts' layout)
        string_cache_path = args.train_string_cache
        if not string_cache_path and args.train_embeddings:
            emb_dir = os.path.dirname(args.train_embeddings)
            candidates = [
                os.path.join(emb_dir, "..", "cache", "string.json"),
                os.path.join(emb_dir, "..", "features", "_string_cache.json"),
            ]
            for c in candidates:
                if os.path.exists(c):
                    string_cache_path = c
                    break
            else:
                print("WARNING: no training STRING cache found (checked: "
                      f"{[os.path.normpath(c) for c in candidates]}). "
                      "PPI features 'n_positive_genes_linking_to_this' and "
                      "'frac_positives_linking' will be 0 for all predictions. "
                      "Pass --train-string-cache to fix this.")

        if string_cache_path and os.path.exists(string_cache_path):
            train_string = json.load(open(string_cache_path))
            print(f"Loaded training STRING cache from {os.path.normpath(string_cache_path)} "
                  f"({len(train_string)} genes)")
        elif string_cache_path:
            print(f"WARNING: --train-string-cache path does not exist: {string_cache_path}. "
                  "PPI positive-linking features will be 0.")

        rows = []
        for i, gene in enumerate(gene_names):
            seq, uid = seqs[i], uid_map[gene]
            dom = lib.domain_features(lib.fetch_interpro(uid))
            go = lib.go_features(lib.fetch_go_text(gene, u))
            partners = lib.fetch_string(gene)
            ppi = lib.ppi_features(gene, partners, known_pos, train_string)
            row = {}
            for k, v in lib.physicochemical_features(seq).items(): row[f"physicochemical__{k}"] = v
            for k, v in lib.sequence_composition_features(seq).items(): row[f"sequence_composition__{k}"] = v
            for k, v in dom.items(): row[f"domain__{k}"] = v
            for k, v in go.items(): row[f"go_term__{k}"] = v
            for k, v in ppi.items(): row[f"ppi_network__{k}"] = v
            rows.append(row)
            if (i+1) % 10 == 0: print(f"    [{i+1:>2}/{N}]")
        bio_df = pd.DataFrame(rows)

        # Align to training column order
        if not args.train_features:
            print("WARNING: --train-features not provided; column order may mismatch.")
            val = pd.concat([plm_df, bio_df], axis=1)
        else:
            train_full = pd.read_csv(args.train_features, nrows=1)
            train_cols = [c for c in train_full.columns if c not in ("Gene","Label","UniProt_ID")]
            combined = pd.concat([plm_df, bio_df], axis=1)
            for c in train_cols:
                if c not in combined.columns: combined[c] = 0
            val = combined[train_cols]

    X = np.nan_to_num(val.values.astype(np.float32))
    print(f"\nFeature matrix: {X.shape}")

    # ── Predict ──
    Xs = pipeline["scaler"].transform(X)
    Xp = pipeline["pca"].transform(Xs)
    prob = pipeline["model"].predict_proba(Xp)[:, 1]
    t = pipeline["optimal_threshold"]
    pred = (prob >= t).astype(int)

    out = pd.DataFrame({
        "Gene": gene_names, "UniProt_ID": [uid_map[g] for g in gene_names],
        "Probability": prob.round(4),
        "Prediction": ["INTERACTOR" if p else "NON-INTERACTOR" for p in pred],
        "Confidence": ["High" if abs(p-t) > 0.3 else "Medium" if abs(p-t) > 0.1 else "Low" for p in prob],
    }).sort_values("Probability", ascending=False)
    out.to_csv(args.output, index=False)

    print(f"\n{'Gene':<15s} {'P(int)':>8s}  {'Prediction':<16s} {'Conf'}")
    print("-" * 50)
    for _, r in out.iterrows():
        print(f"{r['Gene']:<15s} {r['Probability']:>8.4f}  {r['Prediction']:<16s} {r['Confidence']}")
    print(f"\nINTERACTOR: {int(pred.sum())}/{N}, NON-INTERACTOR: {N-int(pred.sum())}/{N}")
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
