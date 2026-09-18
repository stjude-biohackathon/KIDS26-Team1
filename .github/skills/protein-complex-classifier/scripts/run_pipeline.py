#!/usr/bin/env python3
"""
Full pipeline orchestrator: gene lists -> trained classifier.
Runs all steps end-to-end.
"""
import os, sys, argparse, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))


def run(script, extra):
    cmd = ["uv", "run", "python", os.path.join(HERE, script)] + extra
    print(f"\n{'='*70}\n▶ {script} {' '.join(extra)}\n{'='*70}")
    r = subprocess.run(cmd)
    if r.returncode != 0:
        print(f"✗ {script} failed (exit {r.returncode})")
        sys.exit(r.returncode)


def main():
    ap = argparse.ArgumentParser(description="Run the full protein-complex classifier pipeline")
    ap.add_argument("--positive", required=True)
    ap.add_argument("--negative", required=True)
    ap.add_argument("--output-dir", default="output")
    ap.add_argument("--complex-name", default="target complex")
    ap.add_argument("--split", type=float, default=0.3)
    ap.add_argument("--pca-variance", type=float, default=0.95)
    ap.add_argument("--class-weight", default="auto")
    ap.add_argument("--plm-only", action="store_true",
                    help="Skip biological features (PLM embeddings only)")
    ap.add_argument("--skip-embeddings", action="store_true",
                    help="Reuse existing embeddings")
    ap.add_argument("--skip-features", action="store_true",
                    help="Reuse existing biological features")
    args = ap.parse_args()

    print(f"Pipeline for: {args.complex_name}")
    print(f"Positive: {args.positive}\nNegative: {args.negative}\nOutput: {args.output_dir}")

    # Step 1-2: embeddings
    if not args.skip_embeddings:
        run("generate_embeddings.py",
            ["--positive", args.positive, "--negative", args.negative,
             "--output-dir", args.output_dir])

    # Step 3-4: biological features
    if not args.plm_only and not args.skip_features:
        run("extract_biological_features.py", ["--output-dir", args.output_dir])

    # Step 5: train
    train_args = ["--output-dir", args.output_dir, "--split", str(args.split),
                  "--pca-variance", str(args.pca_variance),
                  "--class-weight", args.class_weight, "--best-metric", "F1-Score"]
    if args.plm_only:
        train_args.append("--plm-only")
    run("train_classifier.py", train_args)

    print(f"\n{'='*70}\n✅ PIPELINE COMPLETE\n{'='*70}")
    print(f"  Model:   {args.output_dir}/models/trained_pipeline.pkl")
    print(f"  Results: {args.output_dir}/results/")
    print(f"  Figures: {args.output_dir}/figures/")
    print(f"\n  To predict on new genes:")
    print(f"    uv run python {os.path.join(HERE, 'predict.py')} \\")
    print(f"        --genes <new_gene_list> \\")
    print(f"        --model {args.output_dir}/models/trained_pipeline.pkl \\")
    if not args.plm_only:
        print(f"        --train-features {args.output_dir}/features/full_feature_matrix.csv \\")
        print(f"        --train-embeddings {args.output_dir}/embeddings/combined_feature_matrix.npz \\")
    print(f"        --output predictions.csv")


if __name__ == "__main__":
    main()
