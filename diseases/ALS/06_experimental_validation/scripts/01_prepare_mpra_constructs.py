#!/usr/bin/env python3
"""Build allele-specific MPRA insert sequences from validated Section 4 FASTA records."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
REF_FASTA = ROOT / "diseases/ALS/04_modeling/data/ALS_variants_ref.fa"
ALT_FASTA = ROOT / "diseases/ALS/04_modeling/data/ALS_variants_alt.fa"
OUT = ROOT / "diseases/ALS/06_experimental_validation/results/ALS-S6-R003_MPRA_constructs.tsv"

CANDIDATES = [
    (1, "9:27527365:C:CT", "C", "CT"),
    (2, "1:170074763:G:C", "G", "C"),
    (3, "5:103393614:C:T", "C", "T"),
    (4, "10:2734466:AT:A", "AT", "A"),
]


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def make_insert(sequence: str, allele: str, flank: int = 100) -> str:
    anchor = 1000
    observed = sequence[anchor : anchor + len(allele)]
    if observed != allele:
        raise ValueError(f"allele mismatch at index {anchor}: expected {allele}, observed {observed}")
    return sequence[anchor - flank : anchor] + allele + sequence[anchor + len(allele) : anchor + len(allele) + flank]


def main() -> None:
    ref_records = read_fasta(REF_FASTA)
    alt_records = read_fasta(ALT_FASTA)
    rows: list[dict[str, object]] = []
    for rank, variant, ref, alt in CANDIDATES:
        for label, allele, records in (("REF", ref, ref_records), ("ALT", alt, alt_records)):
            insert = make_insert(records[variant], allele)
            for orientation, sequence in (("forward", insert), ("reverse_complement", reverse_complement(insert))):
                rows.append({
                    "construct_id": f"ALS_S6_{rank:02d}_{label}_{orientation}",
                    "experimental_rank": rank,
                    "variant": variant,
                    "allele_label": label,
                    "allele": allele,
                    "orientation": orientation,
                    "left_flank_bp": 100,
                    "right_flank_bp": 100,
                    "insert_length_bp": len(sequence),
                    "insert_sequence": sequence,
                    "sequence_build": "GRCh38",
                    "source_fasta": str((REF_FASTA if label == "REF" else ALT_FASTA).relative_to(ROOT)),
                })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} constructs to {OUT}")


if __name__ == "__main__":
    main()
