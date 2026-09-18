#!/usr/bin/env python3
"""
Step 3-4: Extract biological features + build full feature matrix.
Requires combined_feature_matrix.npz from generate_embeddings.py.
"""
import os, sys, json, time, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pipeline_lib as lib


def main():
    ap = argparse.ArgumentParser(description="Extract biological features")
    ap.add_argument("--output-dir", default="output", help="Output directory")
    args = ap.parse_args()

    embed_dir = os.path.join(args.output_dir, "embeddings")
    feat_dir = os.path.join(args.output_dir, "features")
    cache_dir = os.path.join(args.output_dir, "cache")
    os.makedirs(feat_dir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)

    # ── Load genes + sequences ──
    npz = np.load(os.path.join(embed_dir, "combined_feature_matrix.npz"), allow_pickle=True)
    gene_names = list(npz["genes"]); y = npz["y"]
    N = len(gene_names)

    with open(os.path.join(cache_dir, "sequences.json")) as f:
        seq_data = json.load(f)
    seqs = [seq_data[g]["sequence"] for g in gene_names]
    uid_map = {g: seq_data[g]["uniprot_id"] for g in gene_names}
    known_positives = set(gene_names[i] for i in range(N) if y[i] == 1)

    print(f"Extracting biological features for {N} proteins...")

    from bioservices import UniProt as UniProtService
    u = UniProtService(verbose=False)

    # ── Caches ──
    def load_cache(name):
        p = os.path.join(cache_dir, name)
        return json.load(open(p)) if os.path.exists(p) else {}

    domain_cache = load_cache("domains.json")
    go_cache = load_cache("go.json")
    string_cache = load_cache("string.json")

    # Fetch missing InterPro/GO/STRING
    print("  Fetching InterPro/GO/STRING (cached where available)...")
    for i, gene in enumerate(gene_names, 1):
        uid = uid_map[gene]
        if gene not in domain_cache:
            domain_cache[gene] = lib.fetch_interpro(uid)
        if gene not in go_cache:
            go_cache[gene] = lib.fetch_go_text(gene, u)
        if gene not in string_cache:
            string_cache[gene] = lib.fetch_string(gene)
        if i % 20 == 0:
            print(f"    [{i:>3}/{N}]")
            json.dump(domain_cache, open(os.path.join(cache_dir, "domains.json"), "w"))
            json.dump(go_cache, open(os.path.join(cache_dir, "go.json"), "w"))
            json.dump(string_cache, open(os.path.join(cache_dir, "string.json"), "w"))
            time.sleep(0.3)
    json.dump(domain_cache, open(os.path.join(cache_dir, "domains.json"), "w"))
    json.dump(go_cache, open(os.path.join(cache_dir, "go.json"), "w"))
    json.dump(string_cache, open(os.path.join(cache_dir, "string.json"), "w"))

    # ── Build per-category feature frames ──
    cats = {"physicochemical": [], "sequence_composition": [],
            "domain": [], "go_term": [], "ppi_network": []}
    for i, gene in enumerate(gene_names):
        seq = seqs[i]
        cats["physicochemical"].append(lib.physicochemical_features(seq))
        cats["sequence_composition"].append(lib.sequence_composition_features(seq))
        cats["domain"].append(lib.domain_features(domain_cache[gene]))
        cats["go_term"].append(lib.go_features(go_cache[gene]))
        cats["ppi_network"].append(lib.ppi_features(gene, string_cache[gene],
                                                    known_positives, string_cache))

    # Save each category + build merged
    bio_parts = []
    for cat, rows in cats.items():
        df = pd.DataFrame(rows)
        df.insert(0, "Gene", gene_names)
        df.to_csv(os.path.join(feat_dir, f"{cat}_features.csv"), index=False)
        print(f"  {cat}: {df.shape[1]-1} features")
        cols = [c for c in df.columns if c != "Gene"]
        bio_parts.append(df[cols].rename(columns={c: f"{cat}__{c}" for c in cols}))

    bio = pd.concat(bio_parts, axis=1)
    bio.insert(0, "Gene", gene_names)
    bio.insert(1, "Label", ["positive" if l == 1 else "negative" for l in y])
    bio.insert(2, "UniProt_ID", [uid_map[g] for g in gene_names])
    bio.to_csv(os.path.join(feat_dir, "all_biological_features.csv"), index=False)
    bio_dim = bio.shape[1] - 3
    print(f"  all_biological: {bio_dim} features")

    # ── Merge with PLM embeddings ──
    plm_df = pd.read_csv(os.path.join(embed_dir, "combined_feature_matrix.csv"))
    plm_only = plm_df.drop(columns=["Gene", "Label", "UniProt_ID"])
    bio_only = bio.drop(columns=["Gene", "Label", "UniProt_ID"])
    full = pd.concat([plm_df[["Gene", "Label", "UniProt_ID"]], plm_only, bio_only], axis=1)
    full.to_csv(os.path.join(feat_dir, "full_feature_matrix.csv"), index=False)

    X = np.nan_to_num(full.drop(columns=["Gene", "Label", "UniProt_ID"]).values.astype(np.float32))
    np.savez(os.path.join(feat_dir, "full_feature_matrix.npz"),
             X=X, y=y, genes=gene_names, plm_dim=plm_only.shape[1], bio_dim=bio_dim)
    print(f"\nfull_feature_matrix: {N} x {X.shape[1]} ({plm_only.shape[1]} PLM + {bio_dim} bio)")
    print(f"Done. Features in {feat_dir}/")


if __name__ == "__main__":
    main()
