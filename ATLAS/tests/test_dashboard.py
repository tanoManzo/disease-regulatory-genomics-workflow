import json
from pathlib import Path
import tempfile
import unittest

from ATLAS import export_dashboard


class ExportTests(unittest.TestCase):
    def make_run(self, root):
        registry = root / "registry.tsv"
        registry.write_text("model_id\tbiosample\tsource\taccessibility_peaks\th3k27ac_peaks\n"
                            "E003\tTest cells\tENCODE\t/private/ENCFF001.bed.gz\t/private/ENCFF002.bed.gz\n")
        training_path = root / "models/E003/training.json"
        training_path.parent.mkdir(parents=True)
        training = {"metrics": {"test": {"roc_auc": 0.9, "average_precision": 0.8}},
                    "counts": {"test": {"0": 2, "1": 1}}, "checkpoint_sha256": "checkpoint"}
        training_path.write_text(json.dumps(training))
        provenance = {"registry": {"path": "/private/registry", "sha256": export_dashboard.sha256_file(registry)},
                      "models": [{"model_id": "E003", "training": {"sha256": export_dashboard.sha256_file(training_path)}}]}
        (root / "genome-task-0000.json").write_text('{"summary":{"prefix":"test"}}')
        (root / "test.run.json").write_text(json.dumps({"provenance": provenance}))
        task = {"chrom": "1", "start": 1, "end": 5}
        plan_path = root / "genome_plan.json"
        plan_path.write_text(json.dumps({"tasks": [task], "provenance_digest": export_dashboard.canonical_digest(provenance)}))
        proof = {"status": "pass", "index": 0, "task": task, "allele_rows": 1,
                 "compressed_bytes": 10, "status_counts": {"ok": 1}, "shard_sha256": {"test.tsv.gz": "abc"}}
        proof["proof_sha256"] = export_dashboard.canonical_digest(proof)
        (root / "final_audit").mkdir()
        (root / "final_audit/task-0000.json").write_text(json.dumps(proof))
        audit = {"status": "pass", "plan_sha256": export_dashboard.sha256_file(plan_path),
                 "validated_chunks": 1, "status_counts": {"ok": 1}, "allele_rows": 1,
                 "compressed_bytes": 10, "selected_models": ["E003"], "scientifically_validated": False,
                 "records_by_chromosome": {"1": 1}, "seconds": 1, "audit_code_sha256": "audit"}
        (root / "genome_final_audit.json").write_text(json.dumps(audit))
        return registry

    def test_proof_checksum_and_status_required(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "proof.json"
            proof = {"status": "pass", "allele_rows": 3}
            proof["proof_sha256"] = export_dashboard.canonical_digest(proof)
            path.write_text(json.dumps(proof))
            self.assertEqual(export_dashboard.checked_proof(path)["allele_rows"], 3)
            proof["allele_rows"] = 4
            path.write_text(json.dumps(proof))
            with self.assertRaises(ValueError):
                export_dashboard.checked_proof(path)

    def test_external_download_urls_are_optional_and_https_only(self):
        self.assertIsNone(export_dashboard.download_url(None, "shard.tsv.gz"))
        self.assertEqual(export_dashboard.download_url("https://example.org/atlas", "shard.tsv.gz"),
                         "https://example.org/atlas/shard.tsv.gz")
        for base in ("javascript:alert(1)", "http://example.org", "https://user:secret@example.org", "https://example.org/?token=x"):
            with self.subTest(base=base), self.assertRaises(ValueError):
                export_dashboard.download_url(base, "file.gz")
        with self.assertRaises(ValueError):
            export_dashboard.download_url("https://example.org", "../private")
        with self.assertRaises(ValueError):
            export_dashboard.download_url(None, "../private")

    def test_successful_public_export_has_no_private_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            registry = self.make_run(root)
            output = root / "site"
            export_dashboard.export(root, registry, output, sample_limit=0)
            snapshot = json.loads((output / "snapshot.json").read_text())
            self.assertEqual(snapshot["rows"], 1)
            self.assertEqual(snapshot["regions"][0]["files"][0]["file"], "test.tsv.gz")
            self.assertFalse(snapshot["scientifically_validated"])
            self.assertEqual(snapshot["models"][0]["dhs"], "ENCFF001")
            for path in output.iterdir():
                self.assertNotIn("/private/", path.read_text())
                self.assertNotIn(str(root), path.read_text())

    def test_training_metadata_changes_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            registry = self.make_run(root)
            with (root / "models/E003/training.json").open("a") as handle:
                handle.write("\n")
            with self.assertRaisesRegex(ValueError, "Training metadata differs"):
                export_dashboard.build_snapshot(root, registry, sample_limit=0)

    def test_registry_changes_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            registry = self.make_run(root)
            registry.write_text("changed")
            with self.assertRaisesRegex(ValueError, "Registry differs"):
                export_dashboard.build_snapshot(root, registry, sample_limit=0)

    def test_mismatched_proof_totals_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            registry = self.make_run(root)
            audit_path = root / "genome_final_audit.json"
            audit = json.loads(audit_path.read_text())
            audit["allele_rows"] = 2
            audit_path.write_text(json.dumps(audit))
            with self.assertRaisesRegex(ValueError, "Proof totals disagree"):
                export_dashboard.build_snapshot(root, registry, sample_limit=0)

    def test_failed_audit_cannot_be_published(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "genome_final_audit.json").write_text('{"status":"failed"}')
            (root / "genome_plan.json").write_text("{}")
            with self.assertRaisesRegex(ValueError, "passed audit"):
                export_dashboard.build_snapshot(root, root / "registry.tsv")


if __name__ == "__main__":
    unittest.main()