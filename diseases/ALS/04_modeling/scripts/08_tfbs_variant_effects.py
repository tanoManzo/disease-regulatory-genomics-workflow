#!/usr/bin/env python3
"""Test whether causal ALS variants alter significant TF-MoDISco/JASPAR sites."""

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


def read_meme(path):
    motifs = {}
    lines = Path(path).read_text().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].startswith("MOTIF "):
            i += 1; continue
        parts = lines[i].split(maxsplit=2)
        motif_id, name = parts[1], parts[2] if len(parts) > 2 else parts[1]
        i += 1
        while i < len(lines) and "letter-probability matrix" not in lines[i]:
            i += 1
        if i >= len(lines): break
        width = int(lines[i].split("w=")[1].split()[0]); i += 1
        matrix = []
        while i < len(lines) and len(matrix) < width:
            if lines[i].strip(): matrix.append([float(x) for x in lines[i].split()[:4]])
            i += 1
        motifs[motif_id] = (name, np.asarray(matrix))
    return motifs


def max_relative(seq, pwm, center=1000):
    eps = 1e-6
    logp = np.log2((pwm + eps) / 0.25)
    rc = logp[::-1][:, [3, 2, 1, 0]]
    minimum = np.minimum(logp, rc).min(axis=1).sum()
    maximum = np.maximum(logp, rc).max(axis=1).sum()
    encode = {"A": 0, "C": 1, "G": 2, "T": 3}
    scores = []
    width = len(pwm)
    for start in range(center - width + 1, center + 1):
        window = seq[start:start + width]
        if len(window) != width or any(base not in encode for base in window):
            continue
        idx = np.array([encode[b] for b in window])
        scores.extend([logp[np.arange(width), idx].sum(), rc[np.arange(width), idx].sum()])
    return (max(scores) - minimum) / (maximum - minimum) if scores and maximum > minimum else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--motif-table", required=True)
    ap.add_argument("--meme", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--ref-fasta", required=True)
    ap.add_argument("--alt-fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rank-out", required=True)
    ap.add_argument("--relative-threshold", type=float, default=0.8)
    args = ap.parse_args()

    table = pd.read_csv(args.motif_table, sep="\t")
    selected = set()
    for _, row in table.iterrows():
        for rank in range(3):
            motif, qval = row.get(f"match{rank}"), row.get(f"qval{rank}")
            if pd.notna(motif) and pd.notna(qval) and float(qval) <= 0.05:
                selected.add(str(motif))
    meme = read_meme(args.meme)
    ref = {r.id: str(r.seq).upper() for r in SeqIO.parse(args.ref_fasta, "fasta")}
    alt = {r.id: str(r.seq).upper() for r in SeqIO.parse(args.alt_fasta, "fasta")}
    variants = pd.read_csv(args.variants, sep="\t")
    variants = variants[variants.name.isin(ref) & variants.name.isin(alt)].copy()
    rows = []
    for v in variants.itertuples(index=False):
        for motif_id in sorted(selected):
            tf, pwm = meme[motif_id]
            rscore, ascore = max_relative(ref[v.name], pwm), max_relative(alt[v.name], pwm)
            rhit, ahit = rscore >= args.relative_threshold, ascore >= args.relative_threshold
            rows.append({"variant": v.name, "motif_id": motif_id, "TF": tf,
                         "causal_regulatory": bool(v.causal_regulatory),
                         "predicted_enhancer_fpr5": bool(v.predicted_enhancer_fpr5),
                         "ref_relative_score": rscore, "alt_relative_score": ascore,
                         "delta_relative_score": ascore-rscore,
                         "variant_in_TFBS": rhit or ahit, "disrupts_TFBS": rhit and not ahit,
                         "creates_TFBS": ahit and not rhit})
    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)
    causal_out = out[out.causal_regulatory]
    rank_rows = []
    for tf, group in out.groupby("TF"):
        cg = group[group.causal_regulatory]
        rank_rows.append({"TF": tf,
            "causal_variants_in_TFBS": cg.loc[cg.variant_in_TFBS, "variant"].nunique(),
            "disrupted_sites": cg.loc[cg.disrupts_TFBS, "variant"].nunique(),
            "created_sites": cg.loc[cg.creates_TFBS, "variant"].nunique(),
            "other_variants_outside_predicted": group.loc[group.variant_in_TFBS & ~group.predicted_enhancer_fpr5, "variant"].nunique(),
            "max_abs_score_change": cg.delta_relative_score.abs().max()})
    rank = pd.DataFrame(rank_rows).sort_values(
        ["causal_variants_in_TFBS", "max_abs_score_change"], ascending=False)
    rank.to_csv(args.rank_out, sep="\t", index=False)
    print(rank.to_string(index=False))
    hit_variants = causal_out.groupby("variant").variant_in_TFBS.any()
    disrupted = causal_out.groupby("variant").disrupts_TFBS.any()
    created = causal_out.groupby("variant").creates_TFBS.any()
    print(f"Causal variants in key TFBS: {hit_variants.sum()}/{len(hit_variants)}")
    print(f"Disrupting: {disrupted.sum()}/{len(disrupted)}; creating: {created.sum()}/{len(created)}")
    per_variant = out.groupby("variant").agg(in_TFBS=("variant_in_TFBS", "any"),
        predicted=("predicted_enhancer_fpr5", "first"), causal=("causal_regulatory", "first"))
    outside = per_variant.in_TFBS & ~per_variant.predicted
    print(f"Other scorable variants in key TFBS outside predicted enhancers: {outside.sum()}/{len(per_variant)}")


if __name__ == "__main__":
    main()
