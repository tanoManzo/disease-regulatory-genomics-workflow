#!/usr/bin/env python3
"""Export ranked TREDNet feature attributions and TF-MoDISco motif matches."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--out-features", required=True)
    ap.add_argument("--out-motifs", required=True)
    args = ap.parse_args()

    run = Path(args.run_dir)
    values = np.load(run / "shap_values_features_TREDNet.npy")
    features = pd.read_csv(args.features, sep="\t")
    if values.ndim == 1:
        values = values[None, :]
    if values.shape[1] != len(features):
        raise ValueError(f"SHAP/features mismatch: {values.shape} versus {len(features)}")
    features["mean_abs_shap"] = np.abs(values).mean(axis=0)
    features["mean_signed_shap"] = values.mean(axis=0)
    features = features.sort_values("mean_abs_shap", ascending=False)
    features.insert(0, "importance_rank", range(1, len(features) + 1))
    features.to_csv(args.out_features, sep="\t", index=False)

    html = run / "tfmodisco_output" / "motifs.html"
    tables = pd.read_html(html)
    motifs = tables[0] if tables else pd.DataFrame()
    motifs.to_csv(args.out_motifs, sep="\t", index=False)
    print(f"Features ranked: {len(features)}; motifs matched: {len(motifs)}")
    print(features.head(20)[["importance_rank", "SRC_TYPE", "TARGET", "CL_NAME",
                              "mean_abs_shap", "mean_signed_shap"]].to_string(index=False))


if __name__ == "__main__":
    main()
