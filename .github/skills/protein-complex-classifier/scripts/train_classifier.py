#!/usr/bin/env python3
"""
Step 5: Train ensemble classifier.
PCA + SMOTE + class weighting + threshold optimization on 9 models.
Requires full_feature_matrix.npz (or combined_feature_matrix.npz with --plm-only).
Saves trained_pipeline.pkl (best model) and all_trained_models.pkl.
"""
import os, sys, argparse, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import joblib

from sklearn.ensemble import (RandomForestClassifier, GradientBoostingClassifier,
    ExtraTreesClassifier, VotingClassifier, StackingClassifier)
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import (confusion_matrix, classification_report, roc_curve, auc,
    precision_recall_curve, average_precision_score, matthews_corrcoef,
    accuracy_score, precision_score, recall_score, f1_score, roc_auc_score)
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from imblearn.over_sampling import SMOTE
from imblearn.ensemble import BalancedBaggingClassifier
warnings.filterwarnings("ignore")

RS = 42


def safe_n_splits(y_arr, requested):
    """Cap StratifiedKFold's n_splits at the smallest class count (and at
    least 2), instead of crashing when a class has fewer members than the
    requested fold count. Prints a note when it reduces the requested value."""
    min_class = int(min(np.bincount(y_arr.astype(int))))
    n = max(2, min(requested, min_class))
    if n < requested:
        print(f"  (note: reduced n_splits {requested}->{n}; smallest class has "
              f"only {min_class} samples)")
    return n


def safe_smote(X_arr, y_arr, random_state, k_requested=5):
    """Run SMOTE if the minority class has enough samples to interpolate
    (>=2); otherwise return the data unchanged rather than crashing. Caps
    k_neighbors at (minority_count - 1), same intent as the original code,
    but never lets it go to 0 or negative."""
    min_class = int(min(np.bincount(y_arr.astype(int))))
    if min_class < 2:
        return X_arr, y_arr
    k = max(1, min(k_requested, min_class - 1))
    sm = SMOTE(random_state=random_state, k_neighbors=k)
    return sm.fit_resample(X_arr, y_arr)


def metrics(yt, yp, yprob):
    tn, fp, fn, tp = confusion_matrix(yt, yp).ravel()
    return {"Accuracy": accuracy_score(yt, yp), "Precision": precision_score(yt, yp, zero_division=0),
            "Recall": recall_score(yt, yp, zero_division=0),
            "Specificity": tn/(tn+fp) if tn+fp else 0, "F1-Score": f1_score(yt, yp, zero_division=0),
            "MCC": matthews_corrcoef(yt, yp), "ROC-AUC": roc_auc_score(yt, yprob),
            "PR-AUC": average_precision_score(yt, yprob), "NPV": tn/(tn+fn) if tn+fn else 0,
            "TP": tp, "FP": fp, "TN": tn, "FN": fn}


def best_threshold(yt, yprob):
    best_f1, best_t = 0, 0.5
    for t in np.arange(0.05, 0.96, 0.01):
        f = f1_score(yt, (yprob >= t).astype(int), zero_division=0)
        if f > best_f1: best_f1, best_t = f, t
    return best_t


def main():
    ap = argparse.ArgumentParser(description="Train ensemble classifier")
    ap.add_argument("--output-dir", default="output")
    ap.add_argument("--plm-only", action="store_true", help="Use combined_feature_matrix (no biological features)")
    ap.add_argument("--split", type=float, default=0.3, help="Test fraction")
    ap.add_argument("--pca-variance", type=float, default=0.95)
    ap.add_argument("--class-weight", default="auto", help="Positive weight or 'auto'")
    ap.add_argument("--no-smote", action="store_true")
    ap.add_argument("--no-threshold-opt", action="store_true")
    ap.add_argument("--best-metric", default="F1-Score")
    args = ap.parse_args()

    feat_dir = os.path.join(args.output_dir, "features")
    embed_dir = os.path.join(args.output_dir, "embeddings")
    model_dir = os.path.join(args.output_dir, "models")
    res_dir = os.path.join(args.output_dir, "results")
    fig_dir = os.path.join(args.output_dir, "figures")
    for d in (model_dir, res_dir, fig_dir): os.makedirs(d, exist_ok=True)

    # ── Load features ──
    if args.plm_only:
        npz = np.load(os.path.join(embed_dir, "combined_feature_matrix.npz"), allow_pickle=True)
    else:
        npz = np.load(os.path.join(feat_dir, "full_feature_matrix.npz"), allow_pickle=True)
    X, y, genes = npz["X"], npz["y"], list(npz["genes"])
    N, D = X.shape
    n_pos, n_neg = int(y.sum()), N - int(y.sum())
    print(f"Dataset: {N} x {D} ({n_pos} pos, {n_neg} neg, ratio 1:{n_neg/n_pos:.1f})")

    weight = round(n_neg / n_pos, 1) if args.class_weight == "auto" else float(args.class_weight)
    cw = {0: 1.0, 1: weight}
    print(f"Class weight: {cw}")

    # ── Split ──
    Xtr, Xte, ytr, yte, itr, ite = train_test_split(
        X, y, np.arange(N), test_size=args.split, stratify=y, random_state=RS)
    genes_te = [genes[i] for i in ite]
    print(f"Train: {len(ytr)} ({int(ytr.sum())} pos), Test: {len(yte)} ({int(yte.sum())} pos)")

    # ── Scale + PCA ──
    scaler = StandardScaler()
    Xtr_s, Xte_s = scaler.fit_transform(Xtr), scaler.transform(Xte)
    pca = PCA(n_components=args.pca_variance, random_state=RS)
    Xtr_p, Xte_p = pca.fit_transform(Xtr_s), pca.transform(Xte_s)
    print(f"PCA: {D} -> {Xtr_p.shape[1]} components ({pca.explained_variance_ratio_.sum()*100:.1f}%)")

    # ── SMOTE ──
    if not args.no_smote:
        Xtr_sm, ytr_sm = safe_smote(Xtr_p, ytr, RS)
        print(f"SMOTE: {int(ytr.sum())} -> {int(ytr_sm.sum())} positive")
    else:
        Xtr_sm, ytr_sm = Xtr_p, ytr

    # ── Classifiers ──
    base = {
        "Random Forest": RandomForestClassifier(n_estimators=500, max_depth=10,
            min_samples_leaf=2, class_weight=cw, random_state=RS, n_jobs=-1),
        "Gradient Boosting": GradientBoostingClassifier(n_estimators=300, max_depth=4,
            learning_rate=0.05, min_samples_leaf=3, subsample=0.8, random_state=RS),
        "SVM (RBF)": SVC(kernel="rbf", C=10.0, gamma="scale", class_weight=cw,
            probability=True, random_state=RS),
        "Logistic Regression": LogisticRegression(C=1.0, class_weight=cw, max_iter=2000,
            solver="lbfgs", random_state=RS),
        "Extra Trees": ExtraTreesClassifier(n_estimators=500, max_depth=10,
            min_samples_leaf=2, class_weight=cw, random_state=RS, n_jobs=-1),
        "MLP": MLPClassifier(hidden_layer_sizes=(256,128), max_iter=1000,
            early_stopping=True, validation_fraction=0.15, alpha=0.001, random_state=RS),
    }
    vote_est = [("rf", base["Random Forest"]), ("svm", base["SVM (RBF)"]),
                ("lr", base["Logistic Regression"]), ("et", base["Extra Trees"]),
                ("mlp", base["MLP"])]
    stacking_cv = safe_n_splits(ytr, 5)
    ensembles = {
        "Soft Voting": VotingClassifier(estimators=vote_est, voting="soft", n_jobs=-1),
        "Stacking (LR meta)": StackingClassifier(estimators=vote_est,
            final_estimator=LogisticRegression(C=1.0, class_weight=cw, max_iter=2000, random_state=RS),
            cv=StratifiedKFold(stacking_cv, shuffle=True, random_state=RS),
            stack_method="predict_proba", n_jobs=-1),
        "Balanced Bagging": BalancedBaggingClassifier(
            estimator=ExtraTreesClassifier(n_estimators=100, max_depth=8,
                class_weight=cw, random_state=RS), n_estimators=20, random_state=RS, n_jobs=-1),
    }
    all_clf = {**base, **ensembles}

    # ── CV for threshold selection ──
    cv_folds = safe_n_splits(ytr, 10)
    print(f"\n{cv_folds}-Fold CV on training set...")
    skf = StratifiedKFold(cv_folds, shuffle=True, random_state=RS)
    thresholds = {}
    for name, clf in all_clf.items():
        yprob_cv = np.zeros(len(ytr))
        for ftr, fval in skf.split(Xtr_p, ytr):
            Xf, yf = (safe_smote(Xtr_p[ftr], ytr[ftr], RS) if not args.no_smote
                      else (Xtr_p[ftr], ytr[ftr]))
            c = clf.__class__(**clf.get_params()) if not isinstance(
                clf, (VotingClassifier, StackingClassifier, BalancedBaggingClassifier)) else clf
            if name in ("Gradient Boosting", "MLP"):
                c.fit(Xf, yf, sample_weight=np.where(yf == 1, weight, 1.0))
            else:
                c.fit(Xf, yf)
            yprob_cv[fval] = c.predict_proba(Xf if False else Xtr_p[fval])[:, 1]
        thresholds[name] = best_threshold(ytr, yprob_cv) if not args.no_threshold_opt else 0.5
        print(f"  {name:<22s} threshold={thresholds[name]:.2f}")

    # ── Retrain on full train, evaluate on test ──
    print("\nEvaluating on test set...")
    test_res = {}
    for name, clf in all_clf.items():
        if name in ("Gradient Boosting", "MLP"):
            clf.fit(Xtr_sm, ytr_sm, sample_weight=np.where(ytr_sm == 1, weight, 1.0))
        else:
            clf.fit(Xtr_sm, ytr_sm)
        yprob = clf.predict_proba(Xte_p)[:, 1]
        t = thresholds[name]
        m = metrics(yte, (yprob >= t).astype(int), yprob)
        test_res[name] = {**m, "threshold": t, "y_prob": yprob,
                          "y_pred": (yprob >= t).astype(int)}
        print(f"  {name:<22s} F1={m['F1-Score']:.3f} Prec={m['Precision']:.3f} "
              f"Rec={m['Recall']:.3f} MCC={m['MCC']:.3f}")

    # ── Best model ──
    best = max(test_res, key=lambda k: test_res[k][args.best_metric])
    print(f"\nBest model: {best} ({args.best_metric}={test_res[best][args.best_metric]:.3f}, "
          f"threshold={test_res[best]['threshold']:.2f})")

    # ── Save pipeline ──
    joblib.dump({"scaler": scaler, "pca": pca, "model": all_clf[best],
                 "optimal_threshold": test_res[best]["threshold"], "best_model_name": best,
                 "class_weight": cw, "feature_dim": D,
                 "plm_only": args.plm_only},
                os.path.join(model_dir, "trained_pipeline.pkl"))
    joblib.dump({"scaler": scaler, "pca": pca, "models": all_clf,
                 "optimal_thresholds": thresholds},
                os.path.join(model_dir, "all_trained_models.pkl"))
    print(f"Saved: {model_dir}/trained_pipeline.pkl + all_trained_models.pkl")

    # ── Results tables ──
    cols = ["threshold","Accuracy","Precision","Recall","Specificity","F1-Score",
            "MCC","ROC-AUC","PR-AUC","NPV","TP","FP","TN","FN"]
    tbl = pd.DataFrame({n: {c: r[c] for c in cols} for n, r in test_res.items()}).T
    tbl.to_csv(os.path.join(res_dir, "test_results.csv"))
    print("\n" + tbl[["threshold","Precision","Recall","F1-Score","MCC","ROC-AUC"]].round(3).to_string())

    # Per-protein test predictions
    pred = pd.DataFrame({"Gene": genes_te,
        "True_Label": ["Positive" if l else "Negative" for l in yte],
        "Predicted": ["Positive" if p else "Negative" for p in test_res[best]["y_pred"]],
        "Probability": test_res[best]["y_prob"].round(4),
        "Correct": yte == test_res[best]["y_pred"]})
    pred.to_csv(os.path.join(res_dir, "test_predictions.csv"), index=False)

    # ── Figures ──
    _figures(test_res, best, yte, ensembles, fig_dir)
    print(f"\nDone. Results in {res_dir}/, figures in {fig_dir}/")


def _figures(test_res, best, yte, ensembles, fig_dir):
    sns.set_theme(style="whitegrid")
    pal = sns.color_palette("husl", len(test_res))
    # Confusion matrix
    fig, ax = plt.subplots(figsize=(6, 5))
    cm = confusion_matrix(yte, test_res[best]["y_pred"])
    sns.heatmap(cm, annot=True, fmt="d", cmap="Greens", ax=ax,
                xticklabels=["Neg","Pos"], yticklabels=["Neg","Pos"], annot_kws={"size":16})
    ax.set_title(f"Confusion Matrix - {best}\n(t={test_res[best]['threshold']:.2f})")
    ax.set_xlabel("Predicted"); ax.set_ylabel("Actual")
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "confusion_matrix.png"), dpi=200)
    # ROC
    fig, ax = plt.subplots(figsize=(8, 6.5))
    for (n, r), c in zip(test_res.items(), pal):
        fpr, tpr, _ = roc_curve(yte, r["y_prob"])
        lw = 2.5 if n in ensembles else 1.5
        ax.plot(fpr, tpr, color=c, lw=lw, label=f"{n} ({auc(fpr,tpr):.3f})")
    ax.plot([0,1],[0,1],"k--", alpha=0.3)
    ax.set_xlabel("FPR"); ax.set_ylabel("TPR"); ax.set_title("ROC Curves - Test Set")
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "roc_curves.png"), dpi=200)
    # PR
    fig, ax = plt.subplots(figsize=(8, 6.5))
    for (n, r), c in zip(test_res.items(), pal):
        p, rc, _ = precision_recall_curve(yte, r["y_prob"])
        ax.plot(rc, p, color=c, lw=2, label=f"{n} (AP={average_precision_score(yte, r['y_prob']):.3f})")
    ax.set_xlabel("Recall"); ax.set_ylabel("Precision"); ax.set_title("PR Curves - Test Set")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "pr_curves.png"), dpi=200)
    plt.close("all")


if __name__ == "__main__":
    main()
