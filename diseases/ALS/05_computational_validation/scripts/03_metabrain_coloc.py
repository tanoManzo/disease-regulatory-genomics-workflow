#!/usr/bin/env python3
"""Candidate-gene coloc sensitivity analysis for ALS GWAS and MetaBrain cortex.

This implements the single-causal-variant approximate Bayes factor model used by
coloc.abf.  It joins the public Project MinE 2021 European ALS GWAS (GRCh37) to
MetaBrain cortex-EUR QTLs (GRCh38) by rsID and allele set, avoiding coordinate-
build ambiguity.  Results are a sensitivity analysis because the one-causal-
variant assumption can be violated at complex loci; the official Project MinE
multi-SNP SMR/HEIDI results are extracted alongside them.
"""

from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path

import pysam
from scipy.special import logsumexp


ROOT = Path(__file__).resolve().parents[4]
RELEASE = ROOT / "data/metabrain/2021-07-23-release"
GWAS = ROOT / "data/als_gwas/project_mine_2021/GWAS/ALS_sumstats_EUR_only.txt.gz"
SMR = ROOT / "data/als_gwas/project_mine_2021/SMR/SMR.multi.ALS.metaBrain.Cortex.summary.mapped.txt"
RESULTS = ROOT / "diseases/ALS/05_computational_validation/results"
LOG = ROOT / "diseases/ALS/05_computational_validation/logs/03_metabrain_coloc.log"
STEM = "2021-07-23-cortex-EUR-80PCs"

# Genes supported or mapped at the two candidate loci with public molecular QTL
# evidence.  C1orf112 is the symbol used by the MetaBrain release (now FIRRM).
LOCUS_GENES = {
    "1:170074763:G:C": ["C1orf112", "KIFAP3", "METTL18", "SCYL3", "SELL"],
    "9:27527365:C:CT": ["C9orf72", "MOB3B", "IFNK"],
}

BASE_COLUMNS = [
    "Gene", "GeneChr", "GenePos", "GeneStrand", "GeneSymbol", "SNP",
    "SNPChr", "SNPPos", "SNPAlleles", "SNPEffectAllele",
    "SNPEffectAlleleFreq", "MetaP", "MetaPN", "MetaPZ", "MetaBeta",
    "MetaSE", "NrDatasets", "DatasetCorrelationCoefficients",
    "DatasetZScores", "DatasetSampleSizes",
]


def complement(allele: str) -> str:
    return allele.translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def compatible_alleles(gwas: set[str], qtl: set[str]) -> bool:
    if gwas == qtl:
        return True
    if all(len(a) == 1 for a in gwas | qtl):
        return {complement(a) for a in gwas} == qtl
    return False


def load_gene_index() -> dict[str, tuple[str, str, int]]:
    wanted = {gene for genes in LOCUS_GENES.values() for gene in genes}
    found: dict[str, tuple[str, str, int]] = {}
    path = RELEASE / f"{STEM}-TopEffects.txt.gz"
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            symbol = row.get("GeneSymbol", "")
            if symbol in wanted:
                found[symbol] = (row["Gene"], row["GeneChr"], int(row["GenePos"]))
    return found


def load_gwas(chromosomes: set[str]) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with gzip.open(GWAS, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["chromosome"] not in chromosomes or not row["rsid"].startswith("rs"):
                continue
            try:
                if float(row["standard_error"]) <= 0:
                    continue
                float(row["beta"])
            except (ValueError, TypeError):
                continue
            old = rows.get(row["rsid"])
            if old is None or float(row["p_value"]) < float(old["p_value"]):
                rows[row["rsid"]] = row
    return rows


def qtl_rows(gene_id: str, chrom: str, pos: int) -> list[dict[str, str]]:
    path = RELEASE / f"{STEM}-chr{chrom}.txt.gz"
    rows: list[dict[str, str]] = []
    with pysam.TabixFile(str(path)) as tabix:
        for line in tabix.fetch(chrom, pos - 1, pos):
            row = dict(zip(BASE_COLUMNS, line.rstrip("\n").split("\t")[: len(BASE_COLUMNS)]))
            if row["Gene"] == gene_id:
                rows.append(row)
    return rows


def log_abf(beta: float, se: float, prior_sd: float) -> float:
    variance = se * se
    prior_variance = prior_sd * prior_sd
    r = prior_variance / (variance + prior_variance)
    z = beta / se
    return 0.5 * (math.log1p(-r) + r * z * z)


def coloc_posteriors(pairs: list[tuple[dict[str, str], dict[str, str]]]) -> dict[str, float]:
    # coloc.abf defaults: p1=p2=1e-4, p12=1e-5; prior SD 0.2 for a
    # case-control log-odds effect and 0.15 for a quantitative-trait effect.
    l1 = [log_abf(float(g["beta"]), float(g["standard_error"]), 0.2) for g, _ in pairs]
    l2 = [log_abf(float(q["MetaBeta"]), float(q["MetaSE"]), 0.15) for _, q in pairs]
    ls1 = float(logsumexp(l1))
    ls2 = float(logsumexp(l2))
    ls12 = float(logsumexp([a + b for a, b in zip(l1, l2)]))
    product_log = ls1 + ls2
    ratio = min(math.exp(ls12 - product_log), 1.0 - 1e-15)
    distinct_log = product_log + math.log1p(-ratio)
    logs = [
        0.0,
        math.log(1e-4) + ls1,
        math.log(1e-4) + ls2,
        math.log(1e-4) + math.log(1e-4) + distinct_log,
        math.log(1e-5) + ls12,
    ]
    denominator = float(logsumexp(logs))
    return {f"pp_h{i}": math.exp(value - denominator) for i, value in enumerate(logs)}


def load_smr() -> dict[str, dict[str, str]]:
    wanted = {gene for genes in LOCUS_GENES.values() for gene in genes}
    rows: dict[str, dict[str, str]] = {}
    with SMR.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Gene_name"] in wanted:
                rows[row["Gene_name"]] = row
    return rows


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    gene_index = load_gene_index()
    gwas = load_gwas({"1", "9"})
    smr = load_smr()
    coloc_output: list[dict[str, str]] = []
    smr_output: list[dict[str, str]] = []

    for candidate, genes in LOCUS_GENES.items():
        for symbol in genes:
            if symbol not in gene_index:
                coloc_output.append({
                    "candidate": candidate, "gene": symbol, "status": "not_in_metabrain_cortex",
                    "qtl_variants": "0", "gwas_rsid_overlap": "0", "allele_matched_variants": "0",
                    **{f"pp_h{i}": "." for i in range(5)}, "best_hypothesis": ".",
                })
                continue
            gene_id, chrom, pos = gene_index[symbol]
            qrows = qtl_rows(gene_id, chrom, pos)
            pairs = []
            overlap = 0
            for qrow in qrows:
                parts = qrow["SNP"].split(":")
                rsid = next((part for part in parts if part.startswith("rs")), "")
                grow = gwas.get(rsid)
                if grow is None:
                    continue
                overlap += 1
                if not compatible_alleles(
                    {grow["effect_allele"].upper(), grow["other_allele"].upper()},
                    {a.upper() for a in qrow["SNPAlleles"].split("/")},
                ):
                    continue
                try:
                    if float(qrow["MetaSE"]) <= 0:
                        continue
                    float(qrow["MetaBeta"])
                except (ValueError, TypeError):
                    continue
                pairs.append((grow, qrow))

            if len(pairs) < 50:
                post = {f"pp_h{i}": float("nan") for i in range(5)}
                status = "insufficient_overlap"
                best = "."
            else:
                post = coloc_posteriors(pairs)
                status = "completed"
                best = max(post, key=post.get).replace("pp_", "").upper()
            coloc_output.append({
                "candidate": candidate,
                "gene": symbol,
                "status": status,
                "qtl_variants": str(len(qrows)),
                "gwas_rsid_overlap": str(overlap),
                "allele_matched_variants": str(len(pairs)),
                **{key: "." if math.isnan(value) else f"{value:.10g}" for key, value in post.items()},
                "best_hypothesis": best,
            })

            if symbol in smr:
                row = smr[symbol]
                smr_output.append({
                    "candidate": candidate,
                    "gene": symbol,
                    "probe_id": row["probeID"],
                    "top_snp_grch38": row["topSNP"],
                    "p_gwas": row["p_GWAS"],
                    "p_eqtl": row["p_eQTL"],
                    "p_smr": row["p_SMR"],
                    "p_smr_multi": row["p_SMR_multi"],
                    "p_heidi": row["p_HEIDI"],
                    "n_snp_heidi": row["nsnp_HEIDI"],
                })

    coloc_path = RESULTS / "ALS-S5-R013_ProjectMinE_MetaBrain_coloc.tsv"
    with coloc_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(coloc_output[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(coloc_output)

    smr_path = RESULTS / "ALS-S5-R014_ProjectMinE_MetaBrain_SMR.tsv"
    with smr_path.open("w", newline="") as handle:
        fields = ["candidate", "gene", "probe_id", "top_snp_grch38", "p_gwas",
                  "p_eqtl", "p_smr", "p_smr_multi", "p_heidi", "n_snp_heidi"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(smr_output)

    completed = sum(row["status"] == "completed" for row in coloc_output)
    LOG.write_text(
        f"gwas_rsids_chr1_chr9={len(gwas)} genes={len(coloc_output)} "
        f"coloc_completed={completed} smr_rows={len(smr_output)}\n"
    )
    print(LOG.read_text().strip())


if __name__ == "__main__":
    main()
