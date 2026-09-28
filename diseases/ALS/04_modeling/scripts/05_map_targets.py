#!/usr/bin/env python3
"""Map prioritized ALS regulatory variants to genes by TSS proximity and MN Hi-C."""

from __future__ import annotations

import argparse
import bisect
import gzip
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


def read_tss(gtf_path: str):
    tss = defaultdict(list)
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            chrom, start, end, strand, attrs = fields[0], int(fields[3]), int(fields[4]), fields[6], fields[8]
            match = re.search(r'gene_name "([^"]+)"', attrs)
            gene = match.group(1) if match else re.search(r'gene_id "([^"]+)"', attrs).group(1)
            pos = start if strand == "+" else end
            tss[chrom].append((pos, gene))
    for chrom in tss:
        tss[chrom].sort()
    return tss


def genes_in(tss, chrom, start, end):
    entries = tss.get(chrom, [])
    positions = [x[0] for x in entries]
    lo, hi = bisect.bisect_left(positions, start + 1), bisect.bisect_right(positions, end)
    return {gene for _, gene in entries[lo:hi]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--loops-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    variants = pd.read_csv(args.variants, sep="\t")
    variants = variants[variants.causal_regulatory.astype(str).str.lower().eq("true")].copy()
    tss = read_tss(args.gtf)

    loops = defaultdict(list)
    for path in sorted(Path(args.loops_dir).glob("*_loops.bedpe.gz")):
        experiment = path.name.split("_")[0]
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) < 6 or f[0] != f[3]:
                    continue
                loops[f[0]].append((int(f[1]), int(f[2]), int(f[4]), int(f[5]), experiment))

    rows = []
    for v in variants.itertuples(index=False):
        entries = tss.get(v.chrom, [])
        positions = [x[0] for x in entries]
        i = bisect.bisect_left(positions, int(v.pos))
        neighbors = entries[max(0, i - 1):i + 1]
        nearest_pos, nearest_gene = min(neighbors, key=lambda x: abs(x[0] - int(v.pos)))
        hic_support = defaultdict(set)
        p0 = int(v.pos) - 1
        for x1, x2, y1, y2, experiment in loops.get(v.chrom, []):
            opposite = None
            if x1 <= p0 < x2:
                opposite = (y1, y2)
            elif y1 <= p0 < y2:
                opposite = (x1, x2)
            if opposite:
                for gene in genes_in(tss, v.chrom, *opposite):
                    hic_support[gene].add(experiment)
        rows.append({
            "priority_rank": v.priority_rank, "variant": v.name, "chrom": v.chrom, "pos": v.pos,
            "nearest_gene": nearest_gene, "nearest_tss_distance_bp": abs(nearest_pos - int(v.pos)),
            "hic_target_genes": ";".join(sorted(hic_support)),
            "hic_supporting_experiments": ";".join(
                f"{gene}:{len(exps)}" for gene, exps in sorted(hic_support.items())),
            "mapping_methods": "proximity;motor_neuron_Hi-C" if hic_support else "proximity",
        })
    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"Mapped {len(out)} causal variants; {(out.hic_target_genes != '').sum()} have a Hi-C target")


if __name__ == "__main__":
    main()
