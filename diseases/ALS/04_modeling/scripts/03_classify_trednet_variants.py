#!/usr/bin/env python3
"""Combine REF/ALT TREDNet scores and classify ALS regulatory variants."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def bool_col(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("true")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref-scores", required=True)
    parser.add_argument("--alt-scores", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--classification", required=True)
    parser.add_argument("--mn-overlap", required=True)
    parser.add_argument("--loci", required=True)
    parser.add_argument("--ld-proxies", required=True)
    parser.add_argument("--focal-tags", required=True)
    parser.add_argument("--model-thresholds", required=True)
    parser.add_argument("--out-variants", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-loci", required=True)
    parser.add_argument("--delta-percentile", type=float, default=95.0)
    args = parser.parse_args()

    ref = pd.read_csv(args.ref_scores, sep="\t")[["name", "trednet_score"]].rename(
        columns={"trednet_score": "ref_score"})
    alt = pd.read_csv(args.alt_scores, sep="\t")[["name", "trednet_score"]].rename(
        columns={"trednet_score": "alt_score"})
    meta = pd.read_csv(args.metadata, sep="\t")
    cls = pd.read_csv(args.classification, sep="\t")
    mn = pd.read_csv(args.mn_overlap, sep="\t")[["name", "MN_H3K27ac", "MN_ATAC", "MN_enhancer"]]
    variants = meta.merge(ref, on="name", how="left").merge(alt, on="name", how="left")
    variants = variants.merge(cls.drop(columns=["chrom", "start", "end", "role"]), on="name", how="left")
    variants = variants.merge(mn, on="name", how="left")

    thresholds: dict[int, float] = {}
    with open(args.model_thresholds) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) == 2:
                thresholds[int(fields[0])] = float(fields[1])
    if 5 not in thresholds:
        raise ValueError("model threshold file does not contain the 5% FPR cutoff")

    variants["region_score"] = variants[["ref_score", "alt_score"]].max(axis=1)
    variants["delta_score"] = variants["alt_score"] - variants["ref_score"]
    variants["abs_delta_score"] = variants["delta_score"].abs()
    scorable = variants["abs_delta_score"].dropna()
    delta_cutoff = float(np.percentile(scorable, args.delta_percentile))
    variants["predicted_enhancer_fpr5"] = variants["region_score"] >= thresholds[5]
    variants["causal_regulatory"] = (
        variants["predicted_enhancer_fpr5"]
        & (variants["abs_delta_score"] >= delta_cutoff)
    )
    variants["direction"] = np.select(
        [variants["delta_score"] > 0, variants["delta_score"] < 0],
        ["gain", "loss"], default="no_change")
    variants["associated_regulatory"] = (
        bool_col(variants["enhancer_ref"]) | bool_col(variants["silencer_ref"])
        | bool_col(variants["MN_enhancer"])
    )
    variants["coding"] = bool_col(variants["coding"])
    proxies = pd.read_csv(args.ld_proxies, sep="\t")
    proxy_to_focal = proxies.groupby("SNP_B")["focal_snp"].agg(lambda x: sorted(set(x)))
    variants["linked_focal_tags"] = variants.apply(
        lambda row: ";".join(proxy_to_focal.get(row["name"], [row["name"]])
                             if row["role"] == "ld_proxy" else [row["name"]]), axis=1)
    variants = variants.sort_values(
        ["causal_regulatory", "abs_delta_score", "region_score"], ascending=[False, False, False])
    variants.insert(0, "priority_rank", range(1, len(variants) + 1))

    n = len(variants)
    summary_rows = [
        ("all_candidates", n, n / n),
        ("scorable", variants.ref_score.notna().sum(), variants.ref_score.notna().mean()),
        ("coding", variants.coding.sum(), variants.coding.mean()),
        ("associated_regulatory", variants.associated_regulatory.sum(), variants.associated_regulatory.mean()),
        ("predicted_enhancer_fpr5", variants.predicted_enhancer_fpr5.sum(), variants.predicted_enhancer_fpr5.mean()),
        ("causal_regulatory", variants.causal_regulatory.sum(), variants.causal_regulatory.mean()),
        ("causal_gain", ((variants.causal_regulatory) & (variants.direction == "gain")).sum(),
         ((variants.causal_regulatory) & (variants.direction == "gain")).mean()),
        ("causal_loss", ((variants.causal_regulatory) & (variants.direction == "loss")).sum(),
         ((variants.causal_regulatory) & (variants.direction == "loss")).mean()),
    ]
    summary = pd.DataFrame(summary_rows, columns=["metric", "n", "fraction_all_candidates"])
    summary["model_fpr5_score_cutoff"] = thresholds[5]
    summary["abs_delta_percentile"] = args.delta_percentile
    summary["abs_delta_cutoff"] = delta_cutoff

    loci = pd.read_csv(args.loci, sep="\t")
    focal = pd.read_csv(args.focal_tags, sep="\t")
    locus_rows = []
    for locus in loci.itertuples(index=False):
        locus_tags = set(focal.loc[("chr" + focal.chrom.astype(str) == locus.chrom)
                                   & (focal.pos - 1 >= locus.start)
                                   & (focal.pos - 1 < locus.end), "snp"])
        within = variants[variants.linked_focal_tags.str.split(";").map(
            lambda tags: bool(locus_tags.intersection(tags)))]
        has_coding = bool(within.coding.any())
        has_causal = bool(within.causal_regulatory.any())
        has_regulatory = bool(within.associated_regulatory.any() | within.predicted_enhancer_fpr5.any())
        if has_coding and has_causal:
            category = "coding + causal regulatory"
        elif has_coding:
            category = "coding only"
        elif has_causal:
            category = "causal regulatory only"
        else:
            category = "other"
        locus_rows.append({"gene": locus.gene, "chrom": locus.chrom, "start": locus.start,
                           "end": locus.end, "n_candidates": len(within),
                           "n_causal_regulatory": int(within.causal_regulatory.sum()),
                           "has_coding": has_coding, "has_any_regulatory": has_regulatory,
                           "category": category,
                           "alternative_cell_type_needed": category == "other"})
    locus_table = pd.DataFrame(locus_rows)
    extra = pd.DataFrame([
        {"metric": "gwas_genes_with_causal_regulatory", "n": int((locus_table.n_causal_regulatory > 0).sum()),
         "fraction_all_candidates": float((locus_table.n_causal_regulatory > 0).mean())},
        {"metric": "gwas_genes_with_no_regulatory", "n": int((~locus_table.has_any_regulatory).sum()),
         "fraction_all_candidates": float((~locus_table.has_any_regulatory).mean())},
    ])
    for col, val in [("model_fpr5_score_cutoff", thresholds[5]),
                     ("abs_delta_percentile", args.delta_percentile),
                     ("abs_delta_cutoff", delta_cutoff)]:
        extra[col] = val
    summary = pd.concat([summary, extra], ignore_index=True)

    for path in [args.out_variants, args.out_summary, args.out_loci]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    variants.to_csv(args.out_variants, sep="\t", index=False)
    summary.to_csv(args.out_summary, sep="\t", index=False)
    locus_table.to_csv(args.out_loci, sep="\t", index=False)
    print(f"5% FPR score cutoff: {thresholds[5]:.6f}")
    print(f"Top {100-args.delta_percentile:.1f}% absolute-delta cutoff: {delta_cutoff:.6f}")
    print(f"Causal regulatory variants: {int(variants.causal_regulatory.sum())}/{n}")


if __name__ == "__main__":
    main()
