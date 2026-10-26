# ATLAS Public Dashboard

The static dashboard lives in [site/index.html](site/index.html). It contains a
timestamped snapshot of the passed audit, not a live Slurm connection. No HPC
credentials, patient data, sample genotypes, local filesystem paths or full score
shards are published. The bounded example contains 500 alleles from the first
exported shard of chr22:20,000,001-25,000,000; it is not a genome-wide search.

## Compare Variant Contexts

The example explorer defaults to **All** contexts, with one row per allele and
model. Search a variant to compare its ten context-specific REF, ALT and delta
scores, or select one context. Sorting and TSV exports use the same filtered
comparison rows. Missing rsIDs are displayed as GRCh38 `chr:position:REF:ALT`
identifiers; the original source ID remains unchanged in exported data.

## Local Preview

Open the HTML directly, or preview the site at `http://127.0.0.1:8765/` while the
local preview server is running. The page has no backend or package installation.
Google Fonts and a pinned Lucide build are external visual dependencies; the
data and controls do not need a remote API. Offline fonts fall back locally.

## Refresh the Snapshot

Run from the repository root on the machine holding the audited run:

```bash
models/TREDNET_v2/.venv/bin/python -m ATLAS.export_dashboard
models/TREDNET_v2/.venv/bin/python -m unittest ATLAS.tests.test_dashboard
```

The exporter requires a passed final audit and all per-chunk proofs, verifies
proof checksums and aggregate counts, matches the frozen plan digest, checks
registry/training metadata hashes, and rehashes the shard supplying the examples.
It does not rerun the full audit or rehash every bulk file. Refresh the audit
first if underlying data changed.

Outputs are `data.js`, `snapshot.json`, `audit-summary.json`, and `shards.tsv`.
Only these small public assets, the site sources, this document, the exporter and
its tests are exempted from the existing ATLAS ignore rule. Models, private run
reports and genome-wide score files remain ignored.

## Publish on GitHub Pages

Review the public snapshot before each deployment. Do not force-add the entire
ATLAS directory, models, or local data.

1. Commit and push the dashboard files and `.github/workflows/atlas-pages.yml`
   to the repository's default branch after review.
2. In repository **Settings > Pages**, choose **GitHub Actions** as the build
   source. Repository permissions and the GitHub plan must permit Pages.
3. Under **Actions**, run **Publish ATLAS dashboard** manually. The workflow
   publishes only `ATLAS/site`, not the repository or the HPC run directory.
4. Use the deployment URL returned by the workflow. For the current repository,
   the expected URL is
   `https://tanomanzo.github.io/disease-regulatory-genomics-workflow/`, subject
   to Pages settings and any custom domain. It becomes available after a
   successful deployment.

This is a repository Pages site, so it occupies the repository's Pages root.
The workflow is manual to avoid publishing future snapshots without review.
Subsequent data refreshes require re-exporting locally, committing the public
snapshot and running the workflow again.

## Publish Full Scores Separately

GitHub Pages is not the host for the 13.49 GB genome-wide atlas. Use an approved
public research data repository or HTTPS object storage, subject to source-data
and model redistribution terms. No public bulk-data hosting is configured yet.
Do not add private or signed URLs, credentials, or expiring access tokens.

After the scored shard files have been uploaded under their original names:

```bash
models/TREDNET_v2/.venv/bin/python -m ATLAS.export_dashboard \
  --download-base https://YOUR-PUBLIC-HOST/atlas-release/
```

The manifest then includes a public URL and audited SHA256 for each shard.
Verify remote availability and checksums before republishing. Region JSON
exports contain every relevant shard filename/hash/URL. Files use 1-based
inclusive anchor intervals and contain REF, ALT and ALT-minus-REF per model.
The site remains explicit that biological validation has not been established.

## Snapshot Evidence

The completed audit covers 619 regions, 8,103 shards and 78,229,218 allele rows
across ten models. There are 78,213,330 scored rows, 15,791 ambiguous contexts and
97 unsupported alleles. The exporter never silently changes those denominators.
Held-out model metrics use test chromosomes 8 and 9. ENCODE ERROR caveats for
HepG2/K562 and the small E021 positive set remain visible.