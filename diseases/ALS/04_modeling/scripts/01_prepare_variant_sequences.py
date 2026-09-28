#!/usr/bin/env python3
"""Prepare variant-centred REF and ALT sequences for ALS TREDNet scoring.

The input is the Section 3 candidate table (GRCh38). Each FASTA record is
exactly 2,001 bp and places the first altered base at zero-based index 1,000.
For indels, the two allele sequences share the same left flank and are cropped
to equal length on the right. Reference alleles are checked against hg38.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pysam


def parse_variant(name: str) -> tuple[str, int, str, str]:
    parts = name.split(":")
    if len(parts) != 4:
        raise ValueError(f"expected CHROM:POS:REF:ALT, got {name!r}")
    chrom, pos, ref, alt = parts
    return chrom.removeprefix("chr"), int(pos), ref.upper(), alt.upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--ref-out", required=True)
    parser.add_argument("--alt-out", required=True)
    parser.add_argument("--metadata-out", required=True)
    parser.add_argument("--focal-tags", help="Section 2 ld_focal.tsv (required for rsID tags)")
    parser.add_argument("--vcf-dir", help="directory containing chrN_EUR/EAS.vcf.gz")
    args = parser.parse_args()

    variants = pd.read_csv(args.variants, sep="\t")
    genome = pysam.FastaFile(args.fasta)
    focal = {}
    if args.focal_tags:
        focal_df = pd.read_csv(args.focal_tags, sep="\t")
        focal = {str(r.snp): (str(r.chrom), int(r.pos)) for r in focal_df.itertuples(index=False)}
        focal.update({f"chr{r.chrom}:{int(r.pos)}": (str(r.chrom), int(r.pos))
                      for r in focal_df.itertuples(index=False)})
    vcfs: dict[tuple[str, str], pysam.VariantFile] = {}

    def resolve(name: str, role: str) -> tuple[str, int, str, str]:
        if name.count(":") == 3:
            return parse_variant(name)
        if role != "tag" or name not in focal or not args.vcf_dir:
            raise ValueError(f"cannot resolve alleles for {name!r}")
        chrom, pos1 = focal[name]
        candidates = []
        for pop in ("EUR", "EAS"):
            path = Path(args.vcf_dir) / f"chr{chrom}_{pop}.vcf.gz"
            if not path.exists():
                continue
            key = (chrom, pop)
            vcfs.setdefault(key, pysam.VariantFile(str(path)))
            for record in vcfs[key].fetch(chrom, pos1 - 1, pos1):
                if record.pos == pos1:
                    for alt in record.alts or ():
                        candidates.append((chrom, pos1, record.ref.upper(), alt.upper(), record.id))
        exact = [x for x in candidates if x[4] == name]
        usable = exact or candidates
        unique = {(c, p, r, a) for c, p, r, a, _ in usable}
        if len(unique) != 1:
            raise ValueError(f"ambiguous/missing 1000G alleles for {name}: {sorted(unique)}")
        return next(iter(unique))
    Path(args.ref_out).parent.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, object]] = []
    with open(args.ref_out, "w") as ref_out, open(args.alt_out, "w") as alt_out:
        for row in variants.itertuples(index=False):
            try:
                chrom, pos1, ref, alt = resolve(str(row.name), str(row.role))
            except ValueError as exc:
                records.append({"name": row.name, "role": row.role, "chrom": row.chrom,
                                "pos": int(row.start) + 1, "ref": "", "alt": "",
                                "status": f"unresolved:{exc}"})
                continue
            fasta_chrom = f"chr{chrom}" if f"chr{chrom}" in genome.references else chrom
            pos0 = pos1 - 1
            observed = genome.fetch(fasta_chrom, pos0, pos0 + len(ref)).upper()
            valid = observed == ref
            status = "ok" if valid else f"reference_mismatch:{observed}"
            if not valid:
                records.append({"name": row.name, "role": row.role, "chrom": fasta_chrom,
                                "pos": pos1, "ref": ref, "alt": alt, "status": status})
                continue

            left = genome.fetch(fasta_chrom, pos0 - 1000, pos0).upper()
            # Fetch enough right-flank sequence for either allele, then crop.
            right = genome.fetch(fasta_chrom, pos0 + len(ref),
                                 pos0 + len(ref) + 2001).upper()
            ref_seq = (left + ref + right)[:2001]
            alt_seq = (left + alt + right)[:2001]
            if len(ref_seq) != 2001 or len(alt_seq) != 2001:
                status = "chromosome_boundary"
            elif set(ref_seq + alt_seq) - set("ACGTN"):
                status = "non_iupac_sequence"
            else:
                ref_out.write(f">{row.name}\n{ref_seq}\n")
                alt_out.write(f">{row.name}\n{alt_seq}\n")
            records.append({"name": row.name, "role": row.role, "chrom": fasta_chrom,
                            "pos": pos1, "ref": ref, "alt": alt, "status": status})

    pd.DataFrame(records).to_csv(args.metadata_out, sep="\t", index=False)
    counts = pd.Series([r["status"] for r in records]).value_counts()
    print(f"Input variants: {len(records)}")
    print(f"Scorable variants: {int(counts.get('ok', 0))}")
    for status, count in counts.items():
        if status != "ok":
            print(f"  {status}: {count}")


if __name__ == "__main__":
    main()
