"""Export a public, path-free dashboard snapshot from a passed atlas audit."""

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
from itertools import islice
import json
from pathlib import Path
from urllib.parse import quote, urlparse

REPO = Path(__file__).resolve().parent.parent


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_proof(path):
    proof = json.loads(path.read_text())
    expected = proof.get("proof_sha256")
    content = {key: value for key, value in proof.items() if key != "proof_sha256"}
    if proof.get("status") != "pass" or canonical_digest(content) != expected:
        raise ValueError(f"Invalid audit proof: {path.name}")
    return proof


def download_url(base, name):
    if Path(name).name != name or "\\" in name:
        raise ValueError("Shard name must not contain a directory")
    if not base:
        return None
    parsed = urlparse(base)
    if parsed.scheme != "https" or not parsed.netloc or parsed.query or parsed.fragment or parsed.username:
        raise ValueError("Download base must be a public HTTPS directory URL")
    return base.rstrip("/") + "/" + quote(name)


def build_snapshot(root, registry, sample_limit=500, download_base=None):
    audit_path = root / "genome_final_audit.json"
    audit = json.loads(audit_path.read_text())
    plan_path = root / "genome_plan.json"
    plan = json.loads(plan_path.read_text())
    if audit.get("status") != "pass" or audit["plan_sha256"] != sha256_file(plan_path):
        raise ValueError("A passed audit matching the frozen plan is required")
    if audit["validated_chunks"] != len(plan["tasks"]):
        raise ValueError("Incomplete audit coverage")
    first_report = json.loads((root / "genome-task-0000.json").read_text())
    frozen = json.loads((root / (first_report["summary"]["prefix"] + ".run.json")).read_text())["provenance"]
    if canonical_digest(frozen) != plan["provenance_digest"]:
        raise ValueError("Frozen provenance does not match the audited plan")
    if sha256_file(registry) != frozen["registry"]["sha256"]:
        raise ValueError("Registry differs from the audited registry")
    if not 0 <= sample_limit <= 2000:
        raise ValueError("Sample limit must be between 0 and 2000")
    download_url(download_base, "check.tsv.gz")
    regions, shards, examples = [], [], []
    statuses = Counter()
    for index, task in enumerate(plan["tasks"]):
        proof = checked_proof(root / "final_audit" / f"task-{index:04d}.json")
        if proof["index"] != index or proof["task"] != task:
            raise ValueError("Proof region mismatch")
        statuses.update(proof["status_counts"])
        regions.append({"index": index, **task, "rows": proof["allele_rows"],
                        "status": proof["status"], "status_counts": proof["status_counts"],
                        "bytes": proof["compressed_bytes"], "shards": len(proof["shard_sha256"]), "files": []})
        for name, digest in proof["shard_sha256"].items():
            url = download_url(download_base, name)
            regions[-1]["files"].append({"file": name, "sha256": digest, "url": url})
            shards.append({"region_index": index, "chrom": task["chrom"], "start": task["start"],
                           "end": task["end"], "file": name, "sha256": digest, "url": url})
            if task["chrom"] == "22" and task["start"] == 20000001 and len(examples) < sample_limit:
                path = root / name
                if sha256_file(path) != digest:
                    raise ValueError("Example shard no longer matches the passed audit")
                with gzip.open(path, "rt") as handle:
                    examples.extend(dict(row) for row in islice(csv.DictReader(handle, delimiter="\t"),
                                                              sample_limit - len(examples)))
    if dict(statuses) != audit["status_counts"] or sum(statuses.values()) != audit["allele_rows"]:
        raise ValueError("Proof totals disagree with final audit")
    if sum(region["bytes"] for region in regions) != audit["compressed_bytes"]:
        raise ValueError("Proof byte totals disagree with final audit")
    models = []
    with registry.open() as handle:
        entries = {row["model_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    for identifier in audit["selected_models"]:
        entry = entries[identifier]
        training_path = root / "models" / identifier / "training.json"
        recorded = next(model for model in frozen["models"] if model["model_id"] == identifier)
        if sha256_file(training_path) != recorded["training"]["sha256"]:
            raise ValueError(f"Training metadata differs from the audited model: {identifier}")
        training = json.loads(training_path.read_text())
        models.append({"id": identifier, "biosample": entry["biosample"], "source": entry["source"],
                       "dhs": Path(entry["accessibility_peaks"]).name.split(".")[0],
                       "h3k27ac": Path(entry["h3k27ac_peaks"]).name.split(".")[0],
                       "roc_auc": training["metrics"]["test"]["roc_auc"],
                       "average_precision": training["metrics"]["test"]["average_precision"],
                       "counts": training["counts"], "checkpoint_sha256": training["checkpoint_sha256"],
                       "audit_warning": identifier in {"E118", "E123"},
                       "small_sample": identifier == "E021"})
    snapshot = {"schema_version": 1,
                "exported_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "audit_as_of": datetime.fromtimestamp(audit_path.stat().st_mtime, timezone.utc).isoformat(timespec="seconds"),
                "assembly": "GRCh38", "audit_status": audit["status"],
                "scientifically_validated": audit["scientifically_validated"],
                "rows": audit["allele_rows"], "chunks": audit["validated_chunks"],
                "bytes": audit["compressed_bytes"], "status_counts": audit["status_counts"],
                "records_by_chromosome": audit["records_by_chromosome"],
                "audit_seconds": audit["seconds"], "plan_sha256": audit["plan_sha256"],
                "audit_code_sha256": audit["audit_code_sha256"], "audit_sha256": sha256_file(audit_path),
                "download_base": download_base, "models": models, "regions": regions,
                "examples": examples, "example_scope": "First exported alleles in chr22:20,000,001-25,000,000",
                "repository": "https://github.com/tanoManzo/disease-regulatory-genomics-workflow"}
    return snapshot, shards


def export(root, registry, output, sample_limit=500, download_base=None):
    snapshot, shards = build_snapshot(root, registry, sample_limit, download_base)
    output.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(snapshot, ensure_ascii=True, separators=(",", ":"), allow_nan=False)
    (output / "data.js").write_text("window.ATLAS_DATA=" + payload + ";\n")
    (output / "snapshot.json").write_text(payload + "\n")
    with (output / "shards.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["region_index", "chrom", "start", "end", "file", "sha256", "url"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(shards)
    summary = {key: value for key, value in snapshot.items() if key not in {"examples", "regions"}}
    (output / "audit-summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(f"Exported {len(snapshot['regions'])} regions, {len(shards)} shards, {len(snapshot['examples'])} example alleles")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=REPO / "ATLAS/runs/pilot")
    parser.add_argument("--registry", type=Path, default=REPO / "ATLAS/pilot.tsv")
    parser.add_argument("--output", type=Path, default=REPO / "ATLAS/site")
    parser.add_argument("--sample-limit", type=int, default=500)
    parser.add_argument("--download-base")
    args = parser.parse_args()
    export(args.root, args.registry, args.output, args.sample_limit, args.download_base)


if __name__ == "__main__":
    main()