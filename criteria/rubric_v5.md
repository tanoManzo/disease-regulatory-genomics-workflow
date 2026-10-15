# Regulatory Genomics Manuscript Evaluation Rubric — v5 

October 8, 2026 · Built on the original two-paper rubric (no caps, no citation lookups)

| Fixed | Value |
| --- | --- |
| Anchors (define 100) | Huang & Ovcharenko 2026 (Sci Adv, composite 95); Li et al. 2025 (Nat Commun, 93) |
| Composite | 0.25·C + 0.20·W + 0.40·N_pts + 0.15·J_pts, rounded |
| N_pts | min(100, N × 100 / 60) |
| J_pts | T1 100 · T2 90 · T3 70 · T4 50 · T5 30 · not publishable 0 |

## Scoring constraints

1. **No caps.** Red flags never change a score. Record each one in the *Red-flag comments* column beside the composite.
2. **No lookups.** Never download or search for a citation, impact factor or journal metric. Use only what the draft prints: author names, title, journal name. Impact factors are stated approximately from the journal name.
3. Cite the draft's own text or figure for every deduction. Score each manuscript independently.

## Criterion 1 — Completeness (0–100%)

Ten components, each scored 0 / 0.25 / 0.5 / 0.75 / 1.0 of its weight.

| # | Component | Weight | 1.0 | 0.5 | 0 |
| --- | --- | --- | --- | --- | --- |
| C1 | Question and hypothesis | 6 | Specific gap, testable hypothesis, disease relevance | Topic stated, gap generic | None |
| C2 | Data assembly and QC | 10 | Named datasets, sample sizes, QC, variant source with thresholds | Datasets named, QC missing | Unspecified |
| C3 | Model development | 12 | Architecture, held-out test, auROC/auPRC, baseline | Model described, one metric, no baseline | None / black box |
| C4 | Variant effect scoring | 10 | Formula, threshold, direction, counts | Scores given, ad hoc thresholds | None |
| C5 | Orthogonal computational validation | 12 | ≥ 3 independent lines with matched controls | 1–2 lines or no controls | None |
| C6 | Experimental / independent-dataset validation | 14 | New wet-lab data or ≥ 2 independent experimental datasets with concordance metrics | One external dataset, or concordance with published data | None |
| C7 | Biological interpretation | 10 | Pathways, cell types, TFs tied to disease; ≥ 2 loci end to end | Enrichment lists, one thin example | None |
| C8 | Candidate nomination and testable follow-up | 8 | Named candidate variants/genes with a stated test | Candidate list, no test | None |
| C9 | Statistics | 8 | Named tests, correction, uncertainty | Tests present, weak correction | Central claims untested |
| C10 | Reporting and reproducibility | 10 | Full methods, limitations, public code and data | Methods complete, code/data not public | Missing sections |

Rules: results must be shown, not promised; placeholders score 0. Credit **concordance with published or external datasets** at 0.5 in C6 even without new experiments. **A draft with C6 = 0 cannot score above 75%** (a floor logic, not a red-flag cap: it is the arithmetic ceiling of the components, stated here for transparency). False statements are not deducted here; they reduce C9/C10 only where they make a component unsupported.

**Calibration bands:** ≥ 80% = all validation layers present with some depth; 70–79% = model, validation and interpretation present, one layer thin; 55–69% = pipeline complete but validation or interpretation shallow; 40–54% = pipeline partly executed or mostly descriptive.

## Criterion 2 — Writing (0–100)

Eight sub-scores 0–4: W1 title/abstract ×3, W2 structure ×4, W3 claim–evidence ×4, W4 precision ×3, W5 statistical reporting ×3, W6 figures/legends ×2, W7 hedging/limitations ×2, W8 mechanics/citations ×2. Score = Σ(sub × factor) / 92 × 100, then deductions:

- −5 generic filler, repetition or list-style Results (once).
- −10 per citation that is implausible **on its face** (author names, title and journal do not fit), max −30.
- −5 per numerical inconsistency inside the draft, max −15.

Bands: 90+ = anchor-grade prose; 80–89 = submission-ready; 70–79 = clear but with imprecise claims or reporting gaps; 60–69 = readable, claim–evidence gaps; < 60 = needs rewriting.

## Criterion 3 — Novelty (0–100%)

Tag every Abstract and Results claim: **N1** novel finding (1.0); **N2** novel method, integration or resource (1.0); **R** replication of known result (0.3); **K** general knowledge (0); **U** unsupported (0). N = 100 × Σw / claims.

Tagging rules used in the calibration:
- A new disease-specific candidate list (variants, genes, loci) is N1 when it is generated and tested in the draft.
- A model re-applied to a new disease with the standard pipeline and no new mechanism is R at best for the method, N1 only for the disease-specific findings.
- Statements of established biology, and results that restate the cited literature, are K or R.
- A claim contradicted by the draft's own figures or tables is U.

## Criterion 4 — Journal tier

Highest tier whose every gate is passed (Completeness, Writing in points; Novelty as the raw N).

| Tier | Typical journals (by name) | Gates |
| --- | --- | --- |
| T1 | Nature Genetics, Cell, Nature, Science | C ≥ 95; N ≥ 70; new large-scale data |
| T2 | Nature Communications, Science Advances, Cell Genomics | C ≥ 85; N ≥ 55; C6 ≥ 0.75; ≥ 3 validation lines; W ≥ 80 |
| T3 | AJHG, Genome Biology, Molecular Psychiatry, NAR | C ≥ 70; **N ≥ 40**; C6 ≥ 0.5; W ≥ 70 |
| T4 | Human Molecular Genetics, Briefings in Bioinformatics, NAR Genomics and Bioinformatics | C ≥ 55; N ≥ 25; C5 ≥ 0.5; W ≥ 60 |
| T5 | Frontiers in Genetics, Scientific Reports, BMC Genomics | C ≥ 40; W ≥ 50 |
| — | Not publishable yet | Fails T5 |

The novelty gate (N ≥ 40 for T3, ≥ 25 for T4) is what separates strong-but-incremental drafts (high C, N in the 30s) from the two T3 drafts.

## Red-flag comments (no effect on score)

| Flag | Comment to record |
| --- | --- |
| Fabricated or implausible citation (face check) | List the citation as printed |
| Statements contradicting the draft's figures/tables | Quote both |
| No held-out test / leakage | Describe |
| Circular validation only | Describe |
| No statistical test for a central claim | Name the claim |
| No code/data statement | Note |
| Apparent overlap with published work | Note, recommend integrity review |

## Interpretation

| Composite | Meaning |
| --- | --- |
| 90–100 | Benchmark-equivalent |
| 75–89 | Close to submission; T3 journal |
| 60–74 | Major revision |
| 40–59 | Substantial work needed |
| < 40 | Early draft |

## Output template

| Criterion | Score | Key evidence |
| --- | --- | --- |
| 1. Completeness | __% | Three lowest components |
| 2. Writing | __/100 | Weakest sub-scores; deductions |
| 3. Novelty | __% | Top N1/N2 claims |
| 4. Journal | Tier — primary / fallback | Blocking gate |
| Red-flag comments | none / list | — |
| **Composite** | **__/100** | Band |
