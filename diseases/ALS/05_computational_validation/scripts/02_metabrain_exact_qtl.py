#!/usr/bin/env python3
"""Query exact ALS candidate variants in the complete local MetaBrain release."""

from __future__ import annotations

import csv
import gzip
import math
import re
from collections import defaultdict
from pathlib import Path

import pysam


ROOT = Path(__file__).resolve().parents[4]
RELEASE = ROOT / "data/metabrain/2021-07-23-release"
RESULTS = ROOT / "diseases/ALS/05_computational_validation/results"
LOG = ROOT / "diseases/ALS/05_computational_validation/logs/02_metabrain_exact_qtl.log"
WINDOW = 1_500_000

CANDIDATES = [
    ("1:170074763:G:C", "1", 170074763, "G", "C", "rs522444"),
    ("5:103393614:C:T", "5", 103393614, "C", "T", "."),
    ("9:27527365:C:CT", "9", 27527365, "C", "CT", "rs11410615"),
    ("10:2734466:AT:A", "10", 2734466, "AT", "A", "."),
]

DATASETS = [
    ("basalganglia", "EUR", 30),
    ("cerebellum", "EUR", 60),
    ("cortex", "AFR", 40),
    ("cortex", "EAS", 30),
    ("cortex", "EUR", 80),
    ("hippocampus", "EUR", 30),
    ("spinalcord", "EUR", 20),
]

BASE_COLUMNS = [
    "Gene", "GeneChr", "GenePos", "GeneStrand", "GeneSymbol", "SNP",
    "SNPChr", "SNPPos", "SNPAlleles", "SNPEffectAllele",
    "SNPEffectAlleleFreq", "MetaP", "MetaPN", "MetaPZ", "MetaBeta",
    "MetaSE", "NrDatasets", "DatasetCorrelationCoefficients",
    "DatasetZScores", "DatasetSampleSizes",
]


def dataset_stem(tissue: str, ancestry: str, pcs: int) -> str:
    return f"2021-07-23-{tissue}-{ancestry}-{pcs}PCs"


def load_thresholds(path: Path) -> dict[str, tuple[float, float]]:
    thresholds: dict[str, tuple[float, float]] = {}
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            thresholds[row["Gene"]] = (
                float(row["PvalueNominalThreshold"]), float(row["qval"])
            )
    return thresholds


def allele_match(observed: str, ref: str, alt: str) -> bool:
    alleles = observed.split("/")
    return len(alleles) == 2 and sorted(alleles) == sorted([ref, alt])


def format_float(value: float | None) -> str:
    if value is None or math.isnan(value):
        return "."
    return f"{value:.12g}"


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    raw_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []

    for tissue, ancestry, pcs in DATASETS:
        stem = dataset_stem(tissue, ancestry, pcs)
        thresholds = load_thresholds(RELEASE / f"{stem}-TopEffects.txt.gz")

        for candidate, chrom, pos, ref, alt, rsid in CANDIDATES:
            data_path = RELEASE / f"{stem}-chr{chrom}.txt.gz"
            tested_rows = 0
            position_rows = 0
            exact_rows = 0
            seen: set[str] = set()

            with pysam.TabixFile(str(data_path)) as tabix:
                start = max(0, pos - WINDOW - 1)
                end = pos + WINDOW
                for line in tabix.fetch(chrom, start, end):
                    tested_rows += 1
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < len(BASE_COLUMNS):
                        raise ValueError(f"Unexpected {len(fields)}-column row in {data_path}")
                    row = dict(zip(BASE_COLUMNS, fields[: len(BASE_COLUMNS)]))
                    if row["SNPChr"] != chrom or int(row["SNPPos"]) != pos:
                        continue
                    position_rows += 1
                    if not allele_match(row["SNPAlleles"], ref, alt):
                        continue
                    if line in seen:
                        continue
                    seen.add(line)
                    exact_rows += 1

                    threshold, gene_qval = thresholds.get(
                        row["Gene"], (float("nan"), float("nan"))
                    )
                    meta_p = float(row["MetaP"])
                    significant = (
                        not math.isnan(threshold)
                        and not math.isnan(gene_qval)
                        and gene_qval <= 0.05
                        and meta_p < threshold
                    )
                    raw_rows.append(
                        {
                            "candidate": candidate,
                            "rsid": rsid,
                            "chrom": chrom,
                            "pos_grch38": str(pos),
                            "ref": ref,
                            "alt": alt,
                            "tissue": tissue,
                            "ancestry": ancestry,
                            "pcs": str(pcs),
                            "gene_id": row["Gene"],
                            "gene_symbol": row["GeneSymbol"],
                            "gene_pos_grch38": row["GenePos"],
                            "gene_strand": row["GeneStrand"],
                            "metabrain_snp": row["SNP"],
                            "metabrain_alleles": row["SNPAlleles"],
                            "effect_allele": row["SNPEffectAllele"],
                            "effect_allele_frequency": row["SNPEffectAlleleFreq"],
                            "meta_p": row["MetaP"],
                            "meta_n": row["MetaPN"],
                            "meta_z": row["MetaPZ"],
                            "meta_beta": row["MetaBeta"],
                            "meta_se": row["MetaSE"],
                            "n_datasets": row["NrDatasets"],
                            "nominal_threshold": format_float(threshold),
                            "gene_qval": format_float(gene_qval),
                            "significant": "yes" if significant else "no",
                            "dataset_correlations": row["DatasetCorrelationCoefficients"],
                            "dataset_zscores": row["DatasetZScores"],
                            "dataset_sample_sizes": row["DatasetSampleSizes"],
                        }
                    )

            audit_rows.append(
                {
                    "candidate": candidate,
                    "tissue": tissue,
                    "ancestry": ancestry,
                    "pcs": str(pcs),
                    "cis_window_rows_scanned": str(tested_rows),
                    "position_match_rows": str(position_rows),
                    "exact_allele_match_rows": str(exact_rows),
                    "source_file": str(data_path.relative_to(ROOT)),
                }
            )

    raw_rows.sort(
        key=lambda r: (
            CANDIDATES.index(next(c for c in CANDIDATES if c[0] == r["candidate"])),
            r["tissue"], r["ancestry"], float(r["meta_p"]), r["gene_symbol"],
        )
    )

    raw_path = RESULTS / "ALS-S5-R010_MetaBrain_exact_QTL.tsv"
    with raw_path.open("w", newline="") as handle:
        fieldnames = list(raw_rows[0].keys()) if raw_rows else [
            "candidate", "rsid", "chrom", "pos_grch38", "ref", "alt",
            "tissue", "ancestry", "pcs", "gene_id", "gene_symbol",
            "gene_pos_grch38", "gene_strand", "metabrain_snp",
            "metabrain_alleles", "effect_allele", "effect_allele_frequency",
            "meta_p", "meta_n", "meta_z", "meta_beta", "meta_se",
            "n_datasets", "nominal_threshold", "gene_qval", "significant",
            "dataset_correlations", "dataset_zscores", "dataset_sample_sizes",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(raw_rows)

    audit_path = RESULTS / "ALS-S5-R011_MetaBrain_query_audit.tsv"
    with audit_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(audit_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)

    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in raw_rows:
        grouped[(row["candidate"], row["tissue"], row["ancestry"])].append(row)
    summary_path = RESULTS / "ALS-S5-R012_MetaBrain_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        fields = [
            "candidate", "tissue", "ancestry", "exact_associations",
            "significant_associations", "genes", "significant_genes", "min_p",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(grouped):
            rows = grouped[key]
            sig = [r for r in rows if r["significant"] == "yes"]
            writer.writerow(
                {
                    "candidate": key[0],
                    "tissue": key[1],
                    "ancestry": key[2],
                    "exact_associations": len(rows),
                    "significant_associations": len(sig),
                    "genes": ";".join(sorted({r["gene_symbol"] for r in rows})),
                    "significant_genes": ";".join(
                        sorted({r["gene_symbol"] for r in sig})
                    ) or ".",
                    "min_p": min(float(r["meta_p"]) for r in rows),
                }
            )

    sig_count = sum(r["significant"] == "yes" for r in raw_rows)
    LOG.write_text(
        f"datasets={len(DATASETS)} candidates={len(CANDIDATES)} "
        f"exact_associations={len(raw_rows)} "
        f"significant_associations={sig_count}\n"
    )
    print(LOG.read_text().strip())


if __name__ == "__main__":
    main()
