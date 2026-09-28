#!/usr/bin/env python3
"""Query exact rsID associations from the public QTLbase2 API."""

from __future__ import annotations

import csv
import json
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
RESULTS = ROOT / "diseases/ALS/05_computational_validation/results"
LOG = ROOT / "diseases/ALS/05_computational_validation/logs/04_qtlbase_query.log"
BASE = "http://www.mulinlab.org/qtlbase"
CANDIDATES = {
    "1:170074763:G:C": "rs522444",
    "9:27527365:C:CT": "rs11410615",
}


def get_json(url: str) -> object:
    request = urllib.request.Request(url, headers={"User-Agent": "ALS-section5-public-query/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def post_json(path: str, values: dict[str, str]) -> object:
    data = urllib.parse.urlencode(values).encode()
    request = urllib.request.Request(
        BASE + path, data=data,
        headers={"User-Agent": "ALS-section5-public-query/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def label(labels: dict, category: str, value: object) -> str:
    entry = labels.get(category, {}).get(str(value), {})
    return entry.get("displayName", str(value)) if isinstance(entry, dict) else str(entry)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    labels = get_json(BASE + "/static/data/sub_table.json")
    rows: list[dict[str, str]] = []
    query_audit: list[dict[str, str]] = []

    for candidate, rsid in CANDIDATES.items():
        initial = post_json("/qtl/query", {"q": rsid.removeprefix("rs"), "t": "0", "d": "0", "g": "hg38"})
        returned = 0
        expected = 0
        for item in initial:
            variant_id = str(item["d"]["id"])
            for cis_trans in ("0", "1"):
                for group in item.get(cis_trans, []):
                    expected += int(group["count"])
                    result = post_json(
                        "/qtl/queryData",
                        {"id": variant_id, "sid": str(group["id"]), "t": "0", "d": "0", "c": cis_trans},
                    )
                    for association in result.get("qtl", []):
                        returned += 1
                        rows.append({
                            "candidate": candidate,
                            "rsid": rsid,
                            "qtlbase_variant_id": variant_id,
                            "association_id": str(association["id"]),
                            "cis_trans": "cis" if cis_trans == "0" else "trans",
                            "qtl_type": label(labels, "xQTL", association["xqtl"]),
                            "tissue": label(labels, "tissue", association["tissue"]),
                            "trait_id": str(association["traitid"]),
                            "effect_allele": str(association.get("effect", "")),
                            "alleles": str(association.get("alt", "")),
                            "effect_size": str(association.get("effectSize", "")),
                            "se": str(association.get("se", "")),
                            "p_value": str(association.get("pval", "")),
                            "fdr": str(association.get("fdr", "")),
                            "dataset_id": str(association.get("metaid", "")),
                        })
        query_audit.append({
            "candidate": candidate,
            "rsid": rsid,
            "status": "completed" if returned == expected else "count_mismatch",
            "expected_associations": str(expected),
            "returned_associations": str(returned),
        })

    trait_meta: dict[str, dict] = {}
    trait_ids = sorted({row["trait_id"] for row in rows}, key=int)
    for start in range(0, len(trait_ids), 100):
        batch = trait_ids[start:start + 100]
        trait_meta.update(post_json("/qtl/qtrait", {"ids": ",".join(batch)}))
    for row in rows:
        meta = trait_meta.get(row["trait_id"], {})
        row["trait_name"] = str(meta.get("hg38_name") or meta.get("name") or row["trait_id"])
        row["trait_symbol"] = str(meta.get("symbol") or "")

    studies = {
        str(study["id"]): study
        for study in post_json("/qtl/getstudies", {})
    }
    for row in rows:
        study = studies.get(row["dataset_id"], {})
        row["dataset_pmid"] = str(study.get("pmid") or "")
        row["dataset_title"] = str(study.get("title") or "")

    raw_path = RESULTS / "ALS-S5-R015_QTLbase_exact_QTL.tsv"
    with raw_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["candidate"], row["qtl_type"], row["cis_trans"])].append(row)
    summary_path = RESULTS / "ALS-S5-R016_QTLbase_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        fields = ["candidate", "qtl_type", "cis_trans", "associations", "tissues", "min_p"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for key in sorted(grouped):
            values = grouped[key]
            writer.writerow({
                "candidate": key[0], "qtl_type": key[1], "cis_trans": key[2],
                "associations": len(values),
                "tissues": ";".join(sorted({row["tissue"] for row in values})),
                "min_p": min(float(row["p_value"]) for row in values),
            })

    audit_path = RESULTS / "ALS-S5-R017_QTLbase_query_audit.tsv"
    with audit_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(query_audit[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(query_audit)

    LOG.write_text(
        f"candidates={len(CANDIDATES)} associations={len(rows)} "
        f"all_counts_match={all(r['status'] == 'completed' for r in query_audit)}\n"
    )
    print(LOG.read_text().strip())


if __name__ == "__main__":
    main()
