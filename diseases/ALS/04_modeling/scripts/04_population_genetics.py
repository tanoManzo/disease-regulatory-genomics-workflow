#!/usr/bin/env python3
"""Annotate TREDNet-prioritized ALS variants with EUR/EAS 1000G frequency."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", required=True)
    ap.add_argument("--vcf-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary-out", required=True)
    args = ap.parse_args()

    variants = pd.read_csv(args.variants, sep="\t")
    handles: dict[tuple[str, str], pysam.VariantFile] = {}

    def af_for(row: pd.Series, pop: str) -> float:
        chrom = str(row.chrom).removeprefix("chr")
        path = Path(args.vcf_dir) / f"chr{chrom}_{pop}.vcf.gz"
        if not path.exists() or not row.ref or not row.alt:
            return np.nan
        key = (chrom, pop)
        handles.setdefault(key, pysam.VariantFile(str(path)))
        for rec in handles[key].fetch(chrom, int(row.pos) - 1, int(row.pos)):
            if rec.pos != int(row.pos) or rec.ref.upper() != str(row.ref).upper():
                continue
            for index, alt in enumerate(rec.alts or ()):
                if alt.upper() == str(row.alt).upper():
                    # Explicit super-population AF; generic AF is not recalculated on subset.
                    af = rec.info.get(f"{pop}_AF")
                    return float(af[index] if isinstance(af, tuple) else af)
        return np.nan

    variants["EUR_AF"] = variants.apply(af_for, axis=1, pop="EUR")
    variants["EAS_AF"] = variants.apply(af_for, axis=1, pop="EAS")
    maximum = variants[["EUR_AF", "EAS_AF"]].max(axis=1, skipna=True)
    variants["frequency_class"] = np.select(
        [maximum >= 0.05, maximum < 0.01], ["common", "rare"], default="low_frequency")
    variants.loc[maximum.isna(), "frequency_class"] = "unavailable"
    variants["population_specific"] = (
        ((variants.EUR_AF >= 0.01) & (variants.EAS_AF < 0.001))
        | ((variants.EAS_AF >= 0.01) & (variants.EUR_AF < 0.001))
    )
    # Ensembl 116 ancestral alleles queried for the causal calls. The chr10
    # repeat deletion cannot be normalized unambiguously and remains unresolved.
    ancestral = {"1:170074763:G:C": "C", "5:103393614:C:T": "C",
                 "9:27527365:C:CT": "C"}
    variants["ancestral_allele"] = variants.name.map(ancestral)
    variants["derived_allele"] = variants.apply(
        lambda r: r.alt if r.ancestral_allele == r.ref else
                  (r.ref if r.ancestral_allele == r.alt else pd.NA), axis=1)
    variants["EUR_DAF"] = variants.apply(
        lambda r: r.EUR_AF if r.derived_allele == r.alt else
                  (1-r.EUR_AF if r.derived_allele == r.ref else np.nan), axis=1)
    variants["EAS_DAF"] = variants.apply(
        lambda r: r.EAS_AF if r.derived_allele == r.alt else
                  (1-r.EAS_AF if r.derived_allele == r.ref else np.nan), axis=1)
    variants["DAF_gt_0_5"] = (variants[["EUR_DAF", "EAS_DAF"]].max(axis=1) > 0.5).astype("boolean")
    variants.loc[variants.derived_allele.isna(), "DAF_gt_0_5"] = pd.NA
    variants.to_csv(args.out, sep="\t", index=False)

    causal = variants[variants.causal_regulatory.astype(str).str.lower().eq("true")]
    rows = []
    for label in ["common", "low_frequency", "rare", "unavailable"]:
        count = int((causal.frequency_class == label).sum())
        rows.append((f"causal_{label}", count, count / len(causal) if len(causal) else np.nan))
    count = int(causal.population_specific.sum())
    rows.append(("causal_population_specific", count, count / len(causal) if len(causal) else np.nan))
    resolved = causal.DAF_gt_0_5.notna()
    rows.append(("causal_DAF_gt_0_5", int(causal.loc[resolved, "DAF_gt_0_5"].sum()),
                 float(causal.loc[resolved, "DAF_gt_0_5"].mean()) if resolved.any() else np.nan))
    summary = pd.DataFrame(rows, columns=["metric", "n", "fraction_causal"])
    summary["limitation"] = "DAF resolved for 3/4 causal variants; chr10 repeat deletion normalization ambiguous"
    summary.to_csv(args.summary_out, sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
