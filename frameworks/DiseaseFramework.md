# Framework of a foundational analysis prefacing a study of regulatory genomics of a particular disease.

## Introductory Query (if using an AI agent)
### Role and Final Goal
You are an autonomous research agent. Your final deliverable is a written scientific manuscript on the regulatory genomics of a specified disease. Every step you take exists to generate the evidence, tables, and figures needed to write that manuscript.

The user will specify the disease, including subtypes if relevant. The user may also provide ICD codes and/or EFO/MONDO terms.

### How to Work
Execute the attached framework in order (Sections 0 through 7). Record Section 0 decisions in the disease workspace's decision log; create a subfolder for each analysis section and the manuscript where you save results and necessary data. Treat it as the analysis plan that produces the manuscript's results. Work through it methodically rather than jumping ahead.

Traceability. Save and log every step: the query run, the database or tool used, the parameters and thresholds chosen, the raw output, and the interpretation. Maintain a running record so any result in the final manuscript can be traced back to its source and reproduced. Assign each result a stable identifier and reference it when it appears in later sections or in the draft.

Disease-agnostic use. The framework applies to any disease. Worked examples, script paths and numbers quoted from completed studies (currently the ALS study in `diseases/ALS_ASTRA/`) illustrate the method only: choose tissues, cell types, cell lines, QTL and TF resources, model tracks and thresholds for the target disease in Sections 0 and 3, adapt the scripts rather than reusing their disease-specific inputs, and never carry findings across diseases.

Precision. Be concise and specific. State exact numbers, thresholds, genome builds, reference panels, tool versions, and cell lines. Do not round away detail or substitute vague summaries for figures.

Handle uncertainty honestly. When the evidence for a step is thin, absent, or conflicting, say so explicitly. State "insufficient information available on [X]" rather than inferring, extrapolating, or filling gaps with plausible-sounding claims. Never fabricate data, citations, effect sizes, or database entries.

Ask when needed. If a step is ambiguous, requires a decision you are not positioned to make (e.g., which cell line to prioritize, which threshold to set), or depends on data you cannot access, pause and ask rather than guessing.

References. Provide a citation for every factual claim, dataset, and prior finding. Use a consistent citation format and compile a complete reference list. Do not cite sources you have not verified.

### Output Requirements for the Manuscript
Write in clean scientific prose suitable for publication.
Remove all markers of AI-generated text. Do not use em-dashes, and avoid formulaic phrasing, filler transitions, and hedging boilerplate.
Structure results to map onto standard manuscript sections (Background/Introduction, Methods, Results, Discussion).
Present quantitative results in tables where the framework calls for them.
Keep claims tied to evidence you generated; flag anything speculative as such.

Before You Begin
Confirm the target disease, then briefly outline your plan for Section 1 and list any databases, tools, or access credentials you will need but do not yet have. Ask any clarifying questions now.

### Local Reference Data (download once, reuse across diseases)

Disease-agnostic reference datasets are already downloaded and verified (core release 2026-09-21; validation additions through 2026-09-25) in `data/` at the project root (`/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/data/`, also reachable via the `/vf/users/...` mount). DO NOT re-download these when the target disease changes; only the disease filter (EFO/MONDO trait term) and disease-specific files change. See `data/README.md` for the full inventory, exact versions, source URLs, and the mapping of each dataset to framework questions. Contents:

- `gwas_catalog/` — full GWAS Catalog dumps: all associations (ontology-annotated), studies, ancestries, EFO trait mappings. Covers Section 2 study/variant/ancestry questions for any disease by filtering on the trait term.
- `efo/` — EFO ontology (.obo/.owl) to map a disease name to trait URIs including child terms.
- `gencode/` — GENCODE v50 GTFs + GRCh38 primary-assembly FASTA (gene bodies, flanking genes, locus lengths).
- `1000genomes/` — Phase 3 GRCh38 phased biallelic SNV+indel VCFs (chr1-22,X) + population panel files, for LD expansion with ancestry-matched populations.
- `vep_cache/` — Ensembl VEP 116 GRCh38 indexed cache tarball (extract next to a VEP install) for coding/noncoding consequence annotation.
- `conservation/` — UCSC hg38 phyloP100way and phastCons100way bigWigs.
- `repeats/` — UCSC hg38 RepeatMasker table (repetitive-element overlap, Section 3).
- `liftover/` — hg19<->hg38 chain files for legacy summary statistics.
- `open_targets/` — Open Targets Platform association datasets + disease/target indexes (parquet), for "linked to other diseases" queries.
- `screen_ccres/` — ENCODE SCREEN cCRE registry V3 (GRCh38, all biosamples): "are GWAS variants in known regulatory elements" (Section 3).
- `roadmap_peaks/` — Roadmap consolidated H3K27ac/H3K27me3/DNase narrowPeaks for ALL 127 epigenomes (hg19; lift to hg38 with liftover/): enhancer/silencer definitions for any disease's tissues (Section 3).
- `ensembl_regulation/` — Ensembl Regulatory Build 116 + motif features (GRCh38).
- `gtex/` — GTEx v10 cis-eQTLs, all 49 tissues: expression-altering variants, tissue selected per disease (Section 3).
- `metabrain/` — complete MetaBrain GRCh38 cis-eQTL release for seven tissue/ancestry groups, with top effects, tabix indexes and the publication supplement containing trans-eQTL, interaction-eQTL and disease enrichment tables. Reuse locally across diseases; redistribution is prohibited by the source terms.
- `fantom5/` — FANTOM5 CAGE enhancer atlas (hg38).
- `jaspar/` — JASPAR 2026 CORE vertebrate TF motifs (Sections 3-4 TFBS analysis).
- `blacklist/` — ENCODE hg38 blacklist v2 exclusion regions.
- `scripts/` — re-runnable download scripts (`wget -c` resumes); `logs/` — download logs.

Disease-specific inputs (e.g., full GWAS summary statistics for a given study, disease epigenomes) belong in `diseases/<disease>/0N_<section>/data/`, not here. Note for downloads on this HPC: compute nodes need the DTN proxy, `export http_proxy=http://dtnNN-e0:3128 https_proxy=$http_proxy` (NN = 17-27 or 30).


## 0. Study question, evidence route, and decision rules

Before collecting results, state a testable question about **how a regulatory element works** in the disease-relevant context. Name the sequence feature or element class, the expected direction of activity, the relevant cell type or condition, and the observation that would refute the hypothesis. Disease-variant prioritization may be an application, but the central contribution should establish a regulatory mechanism, a reusable functional resource, a predictive capability, or a tested generalization. Record the primary outcome, comparison groups, candidate-selection rules, and analysis plan before inspecting the corresponding outcome data; label decisions made afterward as exploratory.

Choose the evidence route supported by the project's resources and keep its requirements visible throughout Sections 1–7:

| Route | Core evidence to generate | Independent test of the claim |
| --- | --- | --- |
| Functional screen | A sufficiently powered, quality-controlled MPRA, STARR-seq, CRISPR perturbation, or comparable element-level dataset | Replication and targeted follow-up in a relevant context |
| Predictive model | Sequence-to-activity or element-function predictions, held-out quantitative evaluation, and interpretable sequence determinants | Experimental testing of predictions or a defensible evolutionary test |
| De novo design | Prospective sequence designs with a prespecified prediction and comparator | Measured activity of designed sequences against matched controls |
| Comparative analysis | A harmonized multi-context or multi-species functional dataset and predictive comparisons | Held-out contexts/species, evolutionary evidence, or direct functional testing |

The route is a planning choice, not a claim of completed work. Public assays assembled at scale can support a resource analysis when their provenance and comparability are demonstrated; reanalysis of one public dataset is not a newly generated screen. A classifier trained on chromatin labels is not automatically a quantitative sequence-to-activity model. If experiments cannot be completed, finish and report the computational study honestly, identify the untested predictions, and narrow the title and conclusions. Independent published experimental datasets can satisfy the validation route when independence and quantitative concordance are demonstrated; do not describe them as new wet-lab data. The evaluation criteria in `criteria/rubric_v5.md` are the current project benchmark, not a substitute for a journal's instructions or evidence standards.

Maintain an evidence ledger for each headline claim: supporting result IDs, independent dataset or assay, biological context, denominator, statistical test, uncertainty, limitations, and the strongest justified wording. Distinguish **association**, **predicted activity**, **measured reporter activity**, **endogenous perturbation**, and **disease causality**. A GWAS association, QTL, motif match, chromatin contact, or model score alone does not establish that a variant is causal. Use "candidate regulatory variant" until functional perturbation and locus-level evidence justify stronger language. Plan at least two end-to-end worked examples when the evidence permits, including a close or matched negative comparator that tests discrimination rather than illustrating only a positive hit.

Use the v5 rubric as a readiness check at each milestone. The component definitions and weights are:

| Component | Weight | Evidence to assemble | Framework sections |
| --- | --- | --- | --- |
| C1. Question and hypothesis | 6 | Specific gap, testable hypothesis, and disease relevance | 0-1 |
| C2. Data assembly and QC | 10 | Named datasets, sample sizes, QC, and variant sources with thresholds | 2-3 |
| C3. Model development | 12 | Architecture, independent held-out test, auROC, auPRC, and baseline | 4 |
| C4. Variant effect scoring | 10 | Explicit formula, threshold, direction, and counts | 4 |
| C5. Orthogonal computational validation | 12 | At least three independent evidence lines with matched controls | 5 |
| C6. Experimental / independent-dataset validation | 14 | New wet-lab data or at least two independent experimental datasets with concordance metrics | 5-6 |
| C7. Biological interpretation | 10 | Disease-linked pathways, cell types, and TFs, with at least two loci worked end to end | 3-6 |
| C8. Candidate nomination and testable follow-up | 8 | Named candidate variants or genes with a stated test | 4, 6 |
| C9. Statistics | 8 | Named tests, correction, and uncertainty for central claims | 2-6 |
| C10. Reporting and reproducibility | 10 | Full methods, limitations, public code and data, and accurate access restrictions | 0-7 |

Record each as **completed, partial, absent, or not applicable to the chosen route**, with links to actual results. Score each component at 0, 0.25, 0.5, 0.75, or 1.0 of its weight according to the rubric; non-applicability does not remove a component or renormalize the weights. Results must be shown, not promised: placeholders and proposed validation receive no validation credit. C8 does require a stated future test, but that proposal is not completed C6 evidence. One external dataset or concordance with published data receives 0.5 in C6; at least two independent experimental datasets, each with documented independence from training, thresholds and each other, denominators, tested nulls, allele alignment and uncertainty, receive full C6 credit under this framework (Section 5). Report the rubric prose's conservative 0.5 reading as a sensitivity analysis beside the primary score. Cross-context or cross-species generalization remains necessary when the scientific claim requires it, but is not a separate v5 completeness component.

Plan the evidence that the highest-scoring completed study needed, so it is generated rather than retrofitted:

- **Two model families (C3).** If any pretrained component has incompletely documented training exposure (chromosomes, donors, sequences or assays), also train an independent sequence-only model family from scratch on the same elements and chromosome split, and repeat every primary test with it. Full C3 credit requires one family whose entire training data are documented. A pretrained published predictor such as AlphaGenome (Section 4) adds a third comparator with gene-expression and splicing predictions, but it was trained on ENCODE, FANTOM5 and GTEx tracks and does not replace the from-scratch family.
- **Sensitivity of the primary test (C9).** Prespecify, or clearly label as post hoc, the power or minimum detectable effect of the primary test, a dependence-aware interval (for example a locus-cluster bootstrap) and protocol-matched controls (for example frequency-oriented matched alleles for direction tests). Nonsignificance is not equivalence; report what effect sizes the study could have missed.
- **Prediction corroboration (C5).** Test whether predicted elements are supported by assays that did not contribute to training labels (for example CAGE-defined enhancers, TF ChIP-seq and conservation) against controls matched on training label and sequence composition, in addition to annotation tests of disease variants.
- **External experimental benchmarks (C6).** Identify at least two independent experimental datasets early, document their independence in a table, and evaluate both model families against them.
- **Public code (C10).** Use and cite the public repositories for the models and attribution pipeline (TREDNet v2: https://github.com/tanoManzo/TREDNET_v2; DeepFootprinting: https://github.com/tanoManzo/FootPrinting; workflow framework: https://github.com/tanoManzo/disease-regulatory-genomics-workflow), verify that each is public before stating so, and deposit the study-specific scripts there with an archived identifier. Public model and attribution code with public primary data but undeposited study scripts supports 0.75; full credit requires the study scripts to be public.

Apply no red-flag score caps or penalties beyond the rubric's stated component and writing rules. V5 nevertheless explicitly specifies a completeness ceiling of 75% when C6 = 0. Preserve that special rule in the audit and flag its inconsistency with the rubric's "no caps" statement and arithmetic explanation: the other component weights sum to 86%, not 75%. Report both the weighted component sum and the completeness used for scoring when this rule applies; do not silently reinterpret the source rubric. The former 65% experimental-validation/generalization cap no longer applies. Infer journal tiers only from completed evidence and the gates in Section 7, not from plans.

## 1. Background

### Epidemiology and heritability

- What are the incidence and prevalence of the disease? Stratify by age, sex, geography (including ancestry/population where data exist), education, marital status, and income.

- Is the disease heritable, and to what degree? What are estimated narrow-sense and broad-sense heritability? Summarize evidence from twin and family studies and from SNP-based heritability estimates. How good are these estimates? Does the genetic risk architecture differ by ancestry?

### Genetics, pathways, and affected tissues

- What are the key genes associated with the disease? Have these genes been linked to any other diseases?

- Are there any known gene-gene or gene-environment interactions relevant to this disease?

- Which primary biological pathways are implicated?

- Which cell types and tissues are affected?

- Which ENCODE/NIH Roadmap cell lines are most relevant to those cell types and tissues?

- Assemble the evidence for candidate cell contexts before any variant overlap is inspected, with a citation for each line: tissue and cell-type enrichment reported by the largest GWAS (for example MAGMA, S-LDSC or expression-based enrichment); pathology, imaging or physiology showing which tissues, regions or cell types are affected; single-cell or sorted-cell studies of disease tissue; and peripheral tissues implicated in the disease. List each candidate context with its evidence and mark it as a primary disease context or a comparator. Here "cell line" means any ENCODE or Roadmap biosample (tissue, primary cell, in vitro differentiated cell or immortalized line); Section 3 turns this list into the modeling contexts.

### Comorbidity and progression

- What diseases are comorbid with the disease in question? (OpenTargets can be used to identify comorbidities.)

- Does the disease affect human development in children and teenagers?

- Does the disease progress over time? If so, describe the stages of progression.

### Diagnosis

- What is the current standard-of-care diagnostic approach?

- Are there emerging molecular or genomic diagnostics in use or development?

### Treatment

- What treatments are currently available, if any?

- For each treatment, what is the mechanism of action? Does it target protein, RNA, or DNA?

- For protein-targeting treatments, what class of protein is targeted (e.g., cell-surface receptor, enzyme)?

- Do any treatments act on gene expression or regulation?

- Are any treatments in development that target enhancers or silencers?

- Categorize each treatment as approved, off-label, or investigational.

- Is treatment response known to be affected by genetic variation?

### Open questions

- What are the major unanswered questions, current research focuses, key challenges, and most promising directions in the field?


## 2. GWAS

- What is the total number of GWAS studies targeting this disease? What is the ethnic/racial breakdown of these GWAS cohorts? Additionally, stratify these GWAS cohorts by sex and age. Are there known, specific differences in GWAS results between different populations?

- What is the total number of GWAS-associated SNPs? What is the number of associated indels and translocations?

- What are the top genes most strongly associated with this disease, by top genes we mean genes that have been reproduced in multiple GWAS studies at significant level? Now redefine top genes as those with largest effect sizes while being statistically significant. What are the effect sizes of these associations? From now on we will refer to top genes as GWAS genes. Create a table listing these genes, their left- and right-flanking genes, and their locus lengths.

- What are the most extensively studied coding mutations associated with this disease? Also provide key genetic associations from studies performing exome-wide association studies / rare variant burden analysis from WES/WGS cohorts.

- What fraction of GWAS genes have only non-coding mutations?

- What fraction of all GWAS variants is noncoding?

- Are GWAS gene loci (spanning the gene body + flanking regions) longer or shorter than the genome-wide average? Do they cluster on specific chromosomes? Are these genes evolutionarily well or poorly conserved? Have these genes been linked to any other diseases?

- Identify which of the listed GWAS studies provide summary statistics. Perform fine mapping (using SuSiE, FINEMAP, or PolyFun) with these summary statistics. Save the output of fine mapping. Tabulate all the statistically significant loci with lead variant, credible set, effect size, and allele frequency. Define locus boundaries by credible-set span and TAD or recombination structure. List the flanking genes and locus length.
Record loci with several plausible genes, and loci with no supported genes.

## 3. Regulatory Landscape

- Are there any genes whose expression is known to be altered in this disease, if so in which tissue/cell types?

- List the studies profiling regulatory elements associated with this disease. Indicate which regulatory elements were identified in each study. What are the targets of these regulatory elements? Are any of these targets also genes whose expression is known to be altered in this disease? Are there any variants known to impact the activity of these regulatory elements?

- Are any of the GWAS variants located within known regulatory elements?

- Identify all ATAC-seq, DNAse-seq, ChIP-seq, and other gene regulatory datasets profiling disease relevant cell types. Both bulk and single-cell/nucleus datasets are acceptable for this task.

- Select the modeling contexts systematically and record the choice, with its rationale and rejected alternatives, in the decision log before any variant overlap or heritability result in those contexts is inspected:
  1. Query the ENCODE portal (GRCh38, released, human) for DNase-seq, ATAC-seq, H3K27ac and H3K27me3 experiments in the candidate organs and cell types from Section 1 (organ and cell slims), save the raw responses, and build a biosample-by-assay matrix with donor counts, biosample class (tissue, primary cell, in vitro differentiated cell, cell line), donor disease status and age, and audit flags. List matching Roadmap epigenomes from the local consolidated peaks.
  2. Keep a primary context only if it has accessibility and H3K27ac for an enhancer model, plus H3K27me3 for a silencer model, preferably from the same donors. Record donor numbers per assay and flag marks that come from fewer donors than accessibility.
  3. Prefer the affected cell type, then the affected tissue, then the closest proxy, and state each compromise (bulk tissue, donor age or disease cohort, iPSC-derived or cancer cell lines, missing cell types or developmental stages).
  4. Add comparator contexts not expected to carry risk (for example blood, liver or muscle for a neural disease, or neural tissue for a non-neural one), as models or as heritability annotations, so that context specificity can be tested.
  5. With several donors, fix consensus rules in advance (for example peak support by at least half the donors, or the tissue-level union of histone peaks when donor-level calls are too few), and hold out donors with a recorded seed for donor-level evaluation, excluding donors seen by any pretrained model.
  6. Published tissue enrichment may guide the choice. If the study's own heritability ranking is used to choose contexts, do not then present enrichment in the chosen contexts as an independent finding; otherwise report the heritability partitioning against comparators as a confirmation.
  7. List the contexts that could not be modeled (cell types, developmental stages, peripheral tissues, missing assays) as limitations.

  (Worked examples: the ALS study chose iPSC-derived motor neurons plus cortex, astrocyte, myoblast and monocyte comparators from a systematic ENCODE inventory, `diseases/ALS_ASTRA/01_background/scripts/encode_inventory.py`, `03_regulatory_landscape/scripts/encode_select.py` and `build_contexts.py`; the fibromyalgia study, `frameworks/Fibromyalgia-Ivan-v2.pdf`, chose adult dorsolateral prefrontal cortex and caudate nucleus from GWAS brain enrichment and imaging evidence and used Roadmap blood, liver and muscle epigenomes as heritability comparators.)

- Using H3K27ac and H3K27me3 as preliminary marks of enhancers and silencers, respectively, count the number of enhancers and silencers within GWAS gene loci. For the previous query, redefine enhancers and silencers as regions wherein DHS (or ATAC-seq) peaks overlap with H3K27ac or H3k27me3 marks. Now count the number of enhancers and silencers using this updated definition.

- Focusing on GWAS SNPs/indels, calculate the percentage that are coding, enhancer-based, or silencer-based using this preliminary definition of regulatory elements. What percentage falls within repetitive elements (using UCSC definition)?

- Rank tissues and cell types and genomic context (enhancers/silencers/promoters) by heritability enrichment (S-LDSC with chromatin or expression annotations) and by single-cell methods (scDRS, CELLECT, or similar), using bulk tissue and single-cell ATAC/RNA atlases.

- Treat accessibility, H3K27ac, and H3K27me3 as candidate-element annotations, not direct measurements of enhancer or silencer function. The reference paper explicitly finds that H3K27me3 alone is insufficient to identify silencers. Record the biosample, assay, peak definition, quality controls, and uncertainty for each annotation; test whether putative silencers repress and whether putative enhancers activate in a functional assay when the study's claim depends on those labels.

- Build a functional-data inventory before modeling: assay type; measured sequence or genomic interval; cell type, condition, and donor; construct design; biological and technical replicates; controls; measurement scale; quality-control and significance criteria; availability and reuse terms. If combining public screens, harmonize genome builds, allele orientation, activity definitions, and batch effects. Preserve tested nulls separately from untested elements. State which data will train models, select thresholds, provide independent benchmarks, and support biological interpretation so the same observations are not counted twice.

- For a new functional screen, specify the target sequence universe, positive and negative controls, barcodes or guide RNAs, replicates, batch randomization, power or precision target, quality-control exclusions, and confirmatory assay before ordering constructs or analyzing outcomes. Include both predicted positives and matched predicted negatives, and enough variants or elements from each proposed regulatory class to test its mechanism.

## 4. Modeling

- Train TREDNet (v2, two-phase transfer learning) enhancer and silencer models for the cell line(s)/cell type(s) associated with the disease. Inputs per model, built from cell-type-matched ENCODE data (hg38): positives = DNase or ATAC peaks overlapping H3K27ac (enhancer model) or H3K27me3 (silencer model), as 1 kb windows centered on the accessibility peak summit. Build controls from the reusable genome-wide pool of merged 1 kb DHS windows in `models/TREDNET_v2/input_training_data/control/allDHS.merge.nonPromoterExon`, which excludes promoters (GENCODE TSS +/- 2 kb) and exons. For each model, additionally remove control windows overlapping the ENCODE blacklist, the respective histone mark (H3K27ac or H3K27me3), or any positive window; restrict both classes to canonical chromosomes. Randomly sample controls with a recorded seed at an exact 1:2 positive-to-control ratio overall, and report positive/control counts for the training, validation, and test partitions. Train on chr1-6, chr10-22, chrX/Y; validate on chr7; test on chr8-9. Report test ROC AUC, PR AUC, and score thresholds at 10/5/3/1% FPR (used later as candidate-element thresholds). Consider modeling multiple cell lines/donors where appropriate and compare per-model performance. (Implementation: `models/TREDNET_v2/make_input_training_data.py` for inputs, using `--control-pool`, `--ratio 2`, and a recorded `--seed`; `TREDNet_v2.py` for training; `TREDNet_v2_inference.py` for scoring. Pre-trained models from completed studies are available under `models/TREDNET_v2/`; see the `04_modeling` section of a completed disease workspace under `diseases/` for a worked example with recorded inputs and performance.)

- Establish independence before interpreting model performance. Check whether both phases of any pretrained model saw test chromosomes, overlapping windows, homologous or near-duplicate sequences, donor-derived measurements, or benchmark assay outcomes. Freeze the test set and thresholds before evaluation; if independence cannot be demonstrated, label the result an internal check and assemble a new external test. Report held-out ROC AUC and PR AUC with positive prevalence, confidence intervals, calibration or activity correlation where applicable, and performance by GC content, repeats, cell type, and element class. Compare against prespecified simple baselines such as mark overlap, motif score, conservation, or a published predictor (for example AlphaGenome, below) under the same test conditions. If the study claims quantitative sequence-to-activity prediction, fit and evaluate that outcome directly rather than treating a chromatin-label classifier score as measured activity.

- Train a second, independent model family when the pretrained TREDNet first stage has undocumented training exposure. Use one-hot DNA only, the same positive and control elements, the same 2,001-bp expansion and the same chromosome split (chr7 validation, chr8-9 test), random reverse-complement augmentation, early stopping on validation auROC, forward/reverse-complement averaged inference and a recorded seed (GPU nondeterminism stated). Evaluate it on the same test windows and on held-out donors or samples when available, including sample-specific positives that overlap no training positive; derive its FPR thresholds from its own test negatives; score the same variants and reporter elements; and repeat every primary test, candidate call and external benchmark. Report cross-family agreement of allelic effects (Spearman correlation) and which candidates both families support. Agreement on the primary conclusion removes dependence on the undocumented exposure; disagreement is a result to report. (Worked example: `diseases/ALS_ASTRA/04_modeling/scripts/independent_models.py` on a GPU with the TREDNet environment, then `independent_model_tests.py`.)

- Score variants with AlphaGenome as a pretrained, multimodal published predictor (Avsec et al., *Nature*, doi:10.1038/s41586-025-10014-0; preprint doi:10.1101/2025.06.25.661532). It reads 1 Mb of GRCh38 sequence and predicts 5,563 human tracks (ENCODE 4,769, FANTOM5 CAGE 546, GTEx RNA-seq 216, 4D Nucleome contact maps 28; listed in `models/AlphaGenome/AlphaGenome_HomoSapiense_biosamples.csv`) spanning RNA-seq, CAGE, PRO-cap, DNase, ATAC, histone and TF ChIP-seq, splicing and chromatin contacts.
  - Use it for what TREDNet cannot provide: predicted expression effects (`RNA_SEQ`, gene-masked log fold change) for genes in the 1-Mb window as predicted target-gene evidence beside colocalization (Section 5), never in place of it; splicing effects for variants near exons when classifying coding and regulatory loci; and context-matched accessibility, histone and TF predictions as a published comparator for allelic effects.
  - Prespecify before scoring: the scorers; the disease-relevant tracks, chosen from the cell types and tissues selected for the disease in Section 3 by searching the biosample table on ontology CURIE, `biosample_name` and `gtex_tissue` (for example a disease-relevant primary or differentiated cell type for chromatin tracks and the matched GTEx tissue for `RNA_SEQ`), recording which track types are missing for that context and treating the closest proxy as a sensitivity analysis; and one per-variant summary statistic (for example the quantile score in the prespecified track, or the maximum absolute quantile across prespecified tracks). Score the risk variants and the same matched controls with identical settings, convert REF/ALT effects to the study's risk-allele sign convention, repeat the primary tests, and report Spearman agreement with TREDNet and the sequence-only family and which candidates all families support.
  - Interpret within the scorer. Calibrated quantiles are not P values, confidence intervals or probabilities; directional scorers give signed quantiles from -1 to 1, and `RNA_SEQ` quantiles are on a different scale from other scorers; `*_ACTIVE` scorers do not give effect direction; rows for different tracks or genes are not independent replicates; and raw scores are comparable only within one scorer. Predictions do not establish element function, target-gene causality or mechanism.
  - Profile the genes and genetic and epigenetic features each variant affects, using both scores: `raw_score` is the predicted effect size in the scorer's units (log fold change for the directional scorers) and `quantile_score` is its signed rank among common variants for that track. Summarize each variant by gene for `RNA_SEQ` (and splicing scorers), by histone mark for `CHIP_HISTONE`, by transcription factor for `CHIP_TF`, and by assay for `DNASE`, `ATAC`, `CAGE` and `PROCAP`. For each group report the number of tracks and biosamples, median raw score, median quantile, fraction of tracks with a positive quantile (direction consistency), tracks with |quantile| >= 0.99 and the strongest track, first in the prespecified disease-context tracks and then across all tracks as exploratory. Call a gene or feature affected only when the quantile is extreme, the raw effect exceeds a prespecified minimum, the direction is consistent across context tracks, and its |median quantile| ranks high among the same scorer and feature in matched control variants scored identically; official quantiles alone are not sufficient, because a single variant yields thousands of rows (the demo variant gave 761 of 6,307 `RNA_SEQ` rows with |quantile| >= 0.99, and gene-level quantiles near -0.98 with median raw log fold changes below 0.01 in magnitude) and because the direction can reverse by context (one gene's median quantile in the demo was +0.95 across all tracks but -0.9995 in a two-track tissue-restricted subset). Interpret the profile as a predicted mechanism: report genes by gene type and distance, compare their predicted direction with colocalization and allele-matched eQTL direction; check whether active marks (H3K27ac, H3K4me1, H3K4me3, H3K9ac) and repressive marks (H3K27me3, H3K9me3) move in the directions implied by the TREDNet enhancer or silencer call; and compare TFs with altered ChIP predictions with DeepFootprinting motifs and allele-specific JASPAR motif scores at the variant, noting that most of the 1,617 TF ChIP tracks come from cell lines (HepG2 539, K562 306, HEK293 106, GM12878 101) rather than disease tissue. Convergent gene, chromatin and TF evidence nominates a testable mechanism; it is not validation. (Implementation: `uv run python summarize_scores.py <cases>/*/scores.tsv --controls <controls>/*/scores.tsv --ontology <prespecified CURIEs> --output-dir <dir>` in `models/AlphaGenome/`, which writes `gene_summary.tsv`, `feature_summary.tsv` and a `summary.json` with input checksums.)
  - Check independence. AlphaGenome was trained on ENCODE, FANTOM5, GTEx and 4D Nucleome tracks, so FANTOM5 CAGE or ENCODE chromatin overlap cannot corroborate its predictions, and agreement with TREDNet may reflect shared ENCODE training data. Check each benchmark against the training track list before counting it as independent, and record which genomic intervals the `fold-1` checkpoint held out according to upstream documentation; do not call a prediction held out unless its interval is documented as such.
  - Implementation: `models/AlphaGenome/` (see its README). Install with `uv sync --locked --python 3.13` (pinned `alphagenome` 0.10.0 and `alphagenome-research` 0.3.0 at commit `3eecf48e`); supply a GRCh38 FASTA with its `.fai` through `ALPHAGENOME_FASTA` or `--fasta`; authenticate to Hugging Face after accepting the model terms (`uv run hf auth login` or `HF_TOKEN`; never place tokens in files, notebooks or command lines); and run inference on an allocated GPU (there is no CPU fallback). First run `uv run python score_variants.py --check-only --output-dir <dir>` (REF and calibration check only, no model load), then score each variant with `uv run python score_variants.py --chromosome <chr> --position <1-based> --ref <REF> --alt <ALT> --scorers RNA_SEQ ATAC DNASE CHIP_HISTONE CHIP_TF --output-dir <new dir>`. Positions are 1-based GRCh38, indels must be anchored, and there is no liftover or normalization. Each run writes `scores.tsv` (all gene and track rows with `raw_score` and `quantile_score`) and `run.json` (checkpoint, scorer definitions, versions, calibration SHA-256, device) and refuses to overwrite existing outputs. The runner loads the checkpoint for every variant; for a variant set, write a batch script in the disease workspace that reuses `prepare_variant` and the same checkpoint, scorer and calibration settings, loads the model once per GPU job and writes one manifest per batch. `uv run python -m unittest -v test_score_variants` mocks the model, so confirm GPU inference with the dummy variant first (a software check, not a disease association; a `fold-1` run on 2026-10-15 returned 9,512 score rows on `cuda:0`). Report the checkpoint, versions and commit, scorer definitions, interval length, track selection, summary statistic and calibration checksum in Materials and Methods, cite the upstream repositories (https://github.com/google-deepmind/alphagenome_research, https://github.com/google-deepmind/alphagenome), and do not redistribute model weights, credentials or variant-level outputs from restricted data.

- Quantify what the primary test could detect. Compute exact power of the prespecified test across a grid of true effect sizes at the observed background rate (variant-wise values are upper bounds under LD clustering), the minimum effect detectable with 80% power at the nominal and multiplicity-corrected alpha, and, for direction tests, the critical fraction implied by the matched null. Pair variant-wise intervals with a locus-cluster bootstrap that resamples loci together with their matched controls, and compare risk-allele direction with controls oriented to the same minor- or major-allele status rather than a 50% null. Label these analyses prespecified or post hoc. (Worked examples: `04_modeling/scripts/h1_protocol_tests.py` and `h1_sensitivity_power.py` in the same workspace.)

- Using the TREDNet models trained in the previous step, re-count GWAS variants within or affecting predicted regulatory elements. A variant lies in a predicted enhancer or silencer only under a stated region-score cutoff; its predicted allelic effect is the prespecified difference between allele scores from identically processed variant-centered 2001 bp windows. Record reference, alternate, risk, and protective alleles and the sign convention so a positive score has one unambiguous interpretation. Select candidate-effect thresholds using matched background variants or a held-out functional set, never by optimizing and reporting on the same assay. Show the full distribution, sensitivity to thresholds, counts with denominators, and separate enhancer and silencer calls. Compare candidate regulatory variants with all GWAS variants, coding variants, and regulatory-annotation overlaps. Report the fraction of GWAS genes with at least one candidate regulatory variant and the fraction with none, with locus and gene denominators. Classify loci as coding plus candidate regulatory, coding only, candidate regulatory only, or unresolved. An unresolved locus may reflect cell-context mismatch, missing elements, LD, or model error; test those possibilities before assigning a single cause. Use "candidate" or "predicted effect," reserving "causal" for evidence that isolates the allele or element from linked alternatives in a relevant biological context.

- Using DeepFootprinting (SHAP DeepExplainer applied to the trained TREDNet models), identify key transcription factors (TFs), histone modifications, and other genomic marks associated with enhancers and silencers: compute per-nucleotide SHAP contribution scores over the positive sequences, extract recurring high-importance subsequences as motifs with TF-MoDISco, and match them against the versioned JASPAR motif set used for the run with Tomtom to name the TFs. Calculate the fraction of candidate regulatory variants that fall within one of these key TF binding sites (TFBSs). What fraction of these variants disrupts the corresponding binding site? Do any create new TFBSs? (Implementation: `models/DeepExplainer_TREDNet/` — configure and submit with `submit2biowulf.sh` (single model) or `submit2biowulf_list.sh` (list of models); outputs SHAP scores, feature-importance summary plot, and a TF-MoDISco motif report `tfmodisco_output/motif.html` per run.)

- As a complementary approach, rank the key TFs by the number of candidate regulatory variants located within their binding sites. Determine whether any of the “other” GWAS variants overlap these binding sites outside of the predicted enhancers and silencers.

- Compare and contrast candidate regulatory variants in enhancers and silencers, as well as their associated TFBSs.

- How many of the candidate regulatory variants have been experimentally profiled? Are there experimental studies of noncoding variants that are not predicted to be functional? Why might our model be missing these variants, and how can they be recovered?

- Study the population genetics of candidate regulatory variants. How many are common, rare, or population-specific? How many are derived (DAF > 0.5)?

- Create THE LIST of rank-ordered, prioritized candidate regulatory variants for experimental follow-up.

- For each nominated variant or gene, give its identifier, rank, locus, alleles, relevant context, supporting and contradictory result IDs, evidence level, and a specific follow-up test. State the predicted outcome, matched negative comparator, and observation that would refute the proposed mechanism. A candidate list without a stated test is incomplete C8 evidence; neither the nomination nor its proposed test establishes causality.

- Map the identified enhancers/silencers to their target genes using proximity, re2g, HiChIP, eQTLs, and/or any other method you believe to be appropriate.

- Conduct enrichment analysis of candidate regulatory variants using all GWAS associated variants as background. Use enrichment analysis tools such as GREAT.

- Test regulatory grammar against measured function, not only motif occurrence. For candidate gain or loss of TF binding, quantify allele-specific motif scores, spacing or orientation where relevant, TF expression or occupancy in the assay context, and activity changes in motif-mutant or sequence-swap constructs when feasible. Separate motif enrichment in called elements from evidence that the motif explains their activity. Report negative motifs, null perturbations, and multiple-testing correction.

- Test generalization with genuinely held-out cell types, disease states, donors, or species when the claim extends beyond the training context. Keep genome-build conversion, orthology, sequence conservation, assay scale, and batch effects explicit. Compare performance and activity direction across contexts; a change in score alone is not evidence of a functional switch. De novo designed sequences require prospective measurements and matched controls before claiming successful design.

## 5. Computational Validation

1. QTLs
   - Brain: MetaBrain, SingleBrain, PsychENCODE.
   - General: xQTL Serve, QTLbase, xQTL Atlas, GTEx.
   - Proteome QTLs: [ProteomeVariation](http://proteomevariation.org/).

- Have any of the candidate regulatory elements been supported by eQTLs, sQTLs, and/or other QTLs? Assess QTL evidence in relevant contexts and preserve effect direction.

- Plan at least three independent, orthogonal computational evidence lines for C5, such as fine-mapped QTL colocalization, chromatin-contact target links, and TF occupancy or independent disease-context expression evidence. State what each line tests, its matched controls, statistical test, uncertainty, and provenance. Three resources that reuse the same samples or signal are not automatically three independent lines. Record sample and study overlap, distinguish corroboration from validation, and do not count one result twice merely because it appears in multiple portals. If fewer lines are available, report the gap and score the evidence actually shown.

- Corroborate predicted elements with assays that did not create the training labels. Within held-out test windows of each model, compare windows above the 5% FPR threshold with the remaining windows for overlap with FANTOM5 CAGE-defined enhancers, TF ChIP-seq peaks from the disease-matched cell types (for example BrainTF for brain disorders, or ReMap or ChIP-Atlas peaks filtered to the relevant cell types for other diseases), and mean phyloP, stratified by training label and GC quintile (Cochran-Mantel-Haenszel odds ratios; ordinary least squares with stratum fixed effects and robust errors for continuous scores) with correction across all tests. This tests the model's predictions, not allelic effects, and gives C5 lines that are independent of the association data. Do not apply these resources to AlphaGenome predictions, which were trained on FANTOM5 CAGE and ENCODE tracks. (Worked example from the ALS study: `diseases/ALS_ASTRA/05_computational_validation/scripts/c5_prediction_corroboration.py`.)

- Assess target genes by colocalization rather than by the strongest linked eQTL. Use exact allele-matched GWAS and eQTL summary statistics in each fine-mapping window, approximate Bayes factors with stated priors and a prior sensitivity, a prespecified PP.H4 reporting threshold, and approximate conditional analysis on the lead eQTL for genes favouring distinct variants. Report shared-variant support (H4) and distinct-variant support (H3) as model-dependent posteriors under single-causal-variant assumptions, never as proven target genes, and report tissue sample sizes when tissues differ in power. Reproduce posterior vectors with an independent numerical implementation. Use QTL data from the tissues or cell types relevant to the disease (item 1 above; GTEx v10 covers 49 tissues locally) and state when the relevant context is missing or underpowered. (Worked example with brain eQTLs from the ALS study: `diseases/ALS_ASTRA/05_computational_validation/scripts/coloc_metabrain.py`.)

2. MPRA datasets
   - [MPRAbase](http://mprabase.ucsf.edu/app/mprabase).

- Identify relevant MPRA datasets. Do any of the candidate regulatory variants appear in regulatory active regions identified by MPRA?

- Benchmark at the **exact tested element and allele**, with compatible genome build, strand, and cell context. Separate exact allele matches from interval overlap, which cannot validate an allelic effect. For tested variants, report assay coverage, usable sample size, effect direction and magnitude, uncertainty, significance threshold, correlation or directional concordance with intervals, and performance for positives and tested nulls. Compare with at least one simple baseline and disclose any overlap with model training or threshold selection. Stratify by enhancer versus silencer, assay context, and prediction strength when sufficiently powered. The reference paper compares predicted allele effects with several independent MPRA datasets and uses nearby variants with similar GWAS evidence to test whether its model discriminates their functional effects; reproduce that logic where suitable data exist, without importing its Alzheimer's-specific findings.

- For C6, seek at least two independent experimental datasets with quantitative concordance metrics when new wet-lab data are unavailable. Document independence from model training, threshold selection, and each other; multiple files, tissues, or portal exports from the same experiment do not by themselves establish independence. Report concordance separately for each dataset with denominators, tested nulls, effect-allele alignment, and uncertainty. One external dataset or concordance with published data supports the rubric's 0.5 level; full credit requires the evidence specified in the C6 row. Keep C5 computational evidence and C6 experimental measurements distinct, even when an experiment is analyzed computationally.

- Summarize the external experimental datasets in a main-text table with one row per dataset: assay and cell type, accession and submitting project, measured unit with tested, active and significant denominators, tested nulls, allele or element alignment, the concordance metric with 95% CI for each model family, and the independence statement (training, thresholds, cell context and other datasets). Evaluate at least one model family that shares no training data or cell context with any benchmark, so concordance cannot reflect exposure. Report allelic correlation, directional agreement in the upper half of predicted effects, discrimination of significant from nonsignificant effects, element ranking against matched controls and quantitative activity correlation separately; each supports a different use of the predictions.

- For QTLs, require exact variant or justified fine-mapped/colocalized signal and a plausible target gene in the relevant tissue. Report effect-allele alignment and whether the direction agrees with the predicted enhancer or silencer effect. A QTL hit, chromatin contact, or proximity link is supporting evidence, not a direct perturbation of the regulatory element. Report the number of candidates queried, covered, significant, null, and unavailable in each independent dataset.


3. Biobanks
   - All of Us.
   - Million Veteran Program (MVP).
   - UK Biobank.
   - Biobank Japan.
   - Estonian Biobank.
   - FinnGen (bottleneck population).
   - Qatar Biobank.
   - Genes & Health (consanguineous population; useful for studying higher frequencies of ultrarare variants).
   - China Kadoorie Biobank.

- Validation in independent cohort using summary statistics, via odds ratio, using NHGRI GWAS catalogue.

4. Other relevant datasets
   - BrainTF: ChIP-seq across cell types from different regions of the brain.
   - Genetic-disease catalogs: HGNC, DisGeNET, OMIM, Open Targets.
   - Variant–disease associations: [HuGE Calculator](https://hugeamp.org/hugecalculator.html), Open Targets.
   - Population genetics: [MAtCH](https://match.ctglab.nl/#/home), [Neale Lab UK Biobank LDSC heritability browser](https://nealelab.github.io/UKBB_ldsc/h2_browser.html).
   - Transcription Factors: Cistrome
   - Transcriptomics: GEO
   - Proteomics: Pride
   - Clinically Relevant Variants: ClinVaR
   - Cancer Specific: TCGA, Metabric, Cosmic

- Are there any overlaps with previously experimentally validated variants?

- Optional: Polygenic risk score using candidate regulatory variants with relevant human cohort (All of Us, UK Biobank, Biobank Japan, Fingen, etc.). Examine possibility and/or necessity of conducting our own GWAS/PheWAS.


### Shared local-data availability (updated 2026-09-25)

- Already local in `data/`: GTEx v10 cis-eQTLs, the complete MetaBrain 2021-07-23 cis-eQTL release and publication supplement, MPRAbase v4.9.3, all compact SingleBrain v2 top-association files, NIAGADS xQTL Atlas track metadata, consolidated BrainTF peak collections, HGNC, Open Targets, and Neale Lab UK Biobank LDSC topline results. Section 4 regulatory resources are also local. Disease-specific archives (e.g., full GWAS summary statistics) are documented in the corresponding disease workspace under `diseases/`.
- MetaBrain has restricted reuse terms: use the local copy for scientific research and education, do not redistribute it, and cite de Klein et al. (Nature Genetics, 2023; DOI 10.1038/s41588-023-01300-6). All 315 official manifest files passed MD5 verification.
- Not mirrored because access or licensing remains restricted: PsychENCODE primary data, OMIM, DisGeNET licensed releases, and participant-level biobank data. Obtain the required approval/account before use.
- xQTL Serve is unavailable and functionally superseded by xQTL Atlas. Record the per-disease public portal query status in the disease workspace.
- SingleBrain full association statistics (165.7 GB) are not mirrored; the complete top-association subset is local. Large NIAGADS/QTLbase bulk collections should be added only when the relevant QTL types and tissues have been selected.
- The authoritative path/access inventory is `data/validation_resources.tsv`; general SHA-256 checksums are in `data/validation_checksums.sha256`. Downloads are reproducible with `data/scripts/dl_validation_resources.sh` and `data/scripts/dl_metabrain.sh`; MetaBrain retains its official MD5 manifest and supplemental SHA-256 file inside its release directory.

## 6. Experimental and Biological Validation

This section can use newly generated experiments or independent published experimental datasets that test the claim. For an external-data route, apply the Section 5 exact-element/allele and independence requirements, document the original experimental design and reuse terms, and report measured concordance rather than a proposed experiment. Do not imply that new constructs, perturbations, or measurements were produced locally. New wet-lab work follows the design and reporting requirements below.

Turn the Section 4 priority list into a testable set, including predicted enhancer and silencer effects, high-confidence and borderline scores, matched negatives, and at least two loci suitable for end-to-end mechanism studies. Record why each element and assay context was selected, including whether it is an independent discovery or a follow-up to known functional evidence.

Choose an assay that tests the stated claim. MPRA or STARR-seq measures reporter activity of an assayed sequence; allele-swap reporters test sequence effects; CRISPRi/CRISPRa or genome editing tests the endogenous locus; TF perturbation and motif editing test a proposed mechanism. An assay in an unrelated cell line has limited power to establish a disease-relevant mechanism. When feasible, pair a broad screen with targeted endogenous or orthogonal validation and test the predicted target gene, not only reporter signal.

Before data collection, specify the construct or guide design, sequence orientation, allele labels, relevant cell type and disease condition, controls, randomization, replicate and donor structure, power or minimum detectable effect, quality-control exclusions, primary outcome, statistical test, correction for multiple comparisons, and success/failure criteria. Include negative controls and a biologically close comparator at the worked loci. Check off-target effects, viability, and batch effects where relevant. Preserve raw counts, processed activity, scripts, and provenance.

After the experiment, report all prespecified tests, including null and contradictory outcomes. Give the tested and passing denominators, effect sizes with uncertainty, replicate agreement, and exact significance criteria. Distinguish a validated **regulatory activity** from a validated **allelic effect**, **target-gene effect**, or **disease mechanism**; each requires its own evidence. If the work remains with collaborators or only a proposal exists, label it planned, do not include it as a completed result, and state how its absence limits the manuscript's claims and rubric score.


## 7. Manuscript Writing and Assembly

### Objective

Write one complete scientific manuscript that converts the evidence generated in Sections 0 through 6 into a coherent study of disease-associated regulatory variation. Call it submission-ready only after the required evidence, reporting, authorship, and human review are complete. The manuscript must not read like six independent reports or a chronological activity log. It must present a focused scientific argument in which each analysis motivates the next, the results support the conclusions, and the limitations are explicit.

Use `reference_papers/sciadv.adz3323.pdf` as the primary model for organization, pacing, quantitative detail, figure integration, and scientific tone. The style reference is:

> Huang D, Ovcharenko I. Silencer variants are key drivers of gene up-regulation in Alzheimer's disease. *Science Advances*. 2026;12:eadz3323. doi:10.1126/sciadv.adz3323.

Follow its observable writing conventions closely, but do not copy or lightly paraphrase its sentences. Do not import its Alzheimer's disease claims, results, citations, terminology, or interpretations into a manuscript about another disease. Cite the reference paper in the manuscript only if it is scientifically relevant, not merely because it was used as a writing model.

The reference paper's observable evidence sequence is: variant universe and cell-context model (Fig. 1); expression and chromatin patterns by regulatory-locus class (Fig. 2); functional distinctions among classes (Fig. 3); cell-state and perturbation responses (Fig. 4); subtype analysis (Fig. 5); independent allele-effect comparisons with reporter assays (Fig. 6); and two sets of detailed locus examples that compare nearby variants and converging evidence (Figs. 7–8). Adapt this sequence to the actual disease and evidence. The paper uses published functional datasets for key comparisons and explicitly calls for further testing of uncharacterized predictions; it is a writing and evidence-flow model, not proof that a new study has met every v5 rubric gate. V5 uses Huang and Ovcharenko 2026 and Li et al. 2025 as its calibration anchors; use the reference paper for style and the current rubric for scoring.

Write in the reference paper's **Science Advances Research Article** layout, not a Nature Article layout: Title; authors; Abstract; INTRODUCTION; RESULTS with declarative subheadings; DISCUSSION; MATERIALS AND METHODS; a Supplementary Materials list; REFERENCES AND NOTES; Acknowledgments. Do not use a "Main" section, Extended Data or superscript citations. Main displays are Fig. 1 to Fig. 8 and Table 1 onward; supplementary items are cited in lowercase as fig. S1 and table S1. The completed ALS study (`diseases/ALS_ASTRA/07_manuscript/`, edition 5: `scripts/manuscript_v5_source.md`, `build_v5_manuscript.py`, `make_v5_figure1.py`, `make_v5_figures.py`, `audit_v5_manuscript.py`) is a worked example of this layout, its build checks and its audit; reuse its structure, not its disease-specific content.

### Required inputs and readiness audit

Before drafting, inventory all outputs from Sections 0 through 6, including result files, tables, figures, statistical outputs, software logs, functional-assay measurements, datasets, and verified literature references. Create an internal result-to-manuscript map containing:

- The stable result identifier.
- The exact quantitative finding and its units.
- The source file and analysis that produced it.
- The relevant sample size, genome build, tissue or cell type, ancestry, thresholds, and statistical test.
- The intended manuscript section, figure, or table.
- The source citation or dataset citation, where applicable.

Check that values reported in the text, figures, tables, legends, and supplementary material agree exactly. Resolve discrepancies from primary outputs rather than choosing the most convenient value.

Complete a manuscript readiness table keyed to rubric components C1–C10. For each component, point to a result, figure, method, and data/code artifact; mark missing work explicitly. Separately assess writing and novelty only after a draft exists. Do not call an unperformed MPRA, planned CRISPR test, or proposed cross-species analysis validation. Report whether the manuscript follows a screen, model-plus-validation, de novo-design, or comparative route, and state any gate that remains unmet. Use the rubric's journal-tier rules only as an internal target, not a prediction of editorial acceptance.

### V5 draft-only evaluation

Keep manuscript preparation and rubric evaluation separate. Preparation requires verified sources as specified below. During evaluation, use only the frozen draft's printed text, figures, tables, author names, citation titles, and journal names. Do not search for or download citations, impact factors, or journal metrics to determine scores. If approximate impact factors are reported, infer them only from the printed journal name and label them approximate. Score each manuscript independently and cite its own text or figure for every deduction. Store the audit outside the submitted manuscript.

1. **Completeness (C):** Sum the ten weighted component scores in Section 0, applying the explicitly stated C6 = 0 rule with the inconsistency note when relevant. Show every component's fractional score, weighted points, evidence, and missing work. The calibration bands are at least 80% (all validation layers present with some depth), 70-79% (one layer thin), 55-69% (validation or interpretation shallow), and 40-54% (partly executed or mostly descriptive). False statements reduce C9/C10 only where they leave a component unsupported; do not invent a separate completeness deduction.
2. **Writing (W):** Score each sub-score from 0 to 4: W1 title/abstract (weight 3), W2 structure (4), W3 claim-evidence (4), W4 precision (3), W5 statistical reporting (3), W6 figures/legends (2), W7 hedging/limitations (2), and W8 mechanics/citations (2). Compute `W = 100 * sum(sub-score * weight) / 92`, then deduct 5 once for generic filler, repetition, or list-style Results; 10 per citation implausible on its face, up to 30; and 5 per numerical inconsistency within the draft, up to 15. Cite the printed evidence for each deduction. Bands are 90+ anchor-grade, 80-89 submission-ready writing, 70-79 clear with gaps, 60-69 readable with claim-evidence gaps, and below 60 requiring rewriting; a writing band alone does not establish overall submission readiness.
3. **Novelty (N):** Create an internal ledger of every Abstract and Results claim, with its location, evidence, and tag: N1 novel finding (1.0), N2 novel method/integration/resource (1.0), R replication (0.3), K general knowledge (0), or U unsupported (0). Compute `N = 100 * sum(claim weights) / number of claims`. A disease-specific candidate list is N1 only when generated and tested in the draft. Reapplying a standard model pipeline is R at best for the method; separately assess the disease-specific findings. Established biology is K or R, and a claim contradicted by the draft's own figures or tables is U. Do not relabel unsupported predictions as discoveries to raise novelty. Use one unit per Abstract sentence and per RESULTS sentence (figure and table legends and headings are not units); protect abbreviations such as "Fig.", "fig.", "figs." and "s.e." from sentence splitting, and record every unit with its line, tag and rationale in a ledger file.
4. **Journal tier (J):** Choose the highest tier for which every gate in the table below is met. C6 and C5 gates refer to fractional component scores, not weighted points; N gates use raw novelty, not normalized novelty points. Name a primary and fallback journal from the rubric and state the blocking gate. These are internal targets, not editorial acceptance forecasts.

| Tier | Typical journals | Required gates | J_pts |
| --- | --- | --- | --- |
| T1 | Nature Genetics, Cell, Nature, Science | C >= 95; N >= 70; new large-scale data | 100 |
| T2 | Nature Communications, Science Advances, Cell Genomics | C >= 85; N >= 55; C6 >= 0.75; at least 3 validation lines; W >= 80 | 90 |
| T3 | AJHG, Genome Biology, Molecular Psychiatry, NAR | C >= 70; N >= 40; C6 >= 0.5; W >= 70 | 70 |
| T4 | Human Molecular Genetics, Briefings in Bioinformatics, NAR Genomics and Bioinformatics | C >= 55; N >= 25; C5 >= 0.5; W >= 60 | 50 |
| T5 | Frontiers in Genetics, Scientific Reports, BMC Genomics | C >= 40; W >= 50 | 30 |
| Not publishable yet | No tier assigned | Fails T5 | 0 |

Compute `N_pts = min(100, N * 100 / 60)` and `Composite = round(0.25 * C + 0.20 * W + 0.40 * N_pts + 0.15 * J_pts)`. Interpret 90-100 as benchmark-equivalent, 75-89 as close to submission (the rubric labels this band T3), 60-74 as major revision, 40-59 as substantial work needed, and below 40 as early draft. A composite-band label does not override the independently checked journal-tier gates.

Record red flags alongside the composite without an additional score change: fabricated or face-implausible citations (as printed), claims contradicting figures/tables (quote both), no held-out test or leakage, circular validation, central claims without statistical tests, no code/data statement, and apparent overlap with published work requiring integrity review. Apply only the scoring consequences explicitly specified above, not extra red-flag penalties.

Use the rubric's output template for the audit:

| Criterion | Score | Key evidence |
| --- | --- | --- |
| 1. Completeness | __% | Three lowest components |
| 2. Writing | __/100 | Weakest sub-scores; deductions |
| 3. Novelty | __% | Top N1/N2 claims |
| 4. Journal | Tier: primary / fallback | Blocking gate |
| Red-flag comments | none / list | Printed evidence; rubric inconsistency if applicable |
| **Composite** | **__/100** | Band |

Do not delay completion solely because an analysis or datum is unavailable. Produce a complete manuscript using the evidence that exists and identify every material gap explicitly. Use direct statements such as, "No reliable ancestry-stratified estimate was available," followed by the consequence for interpretation. Do not invent a value, conceal the gap, or leave unexplained `TODO`, `TBD`, or empty sections. If a missing analysis is essential to the central conclusion, state that the conclusion remains unestablished and specify the analysis needed to establish it.

### Central claim and narrative

Before writing prose, formulate one primary claim that is directly supported by the strongest results. Add two to four secondary claims that explain the genetic architecture, regulatory mechanisms, model performance, biological interpretation, or validation evidence. Every Results subsection must support one of these claims. Exclude technically correct analyses that do not advance the central narrative, or move them to the supplementary portion of the same file.

Organize the Results by scientific findings rather than by framework section number. A typical narrative is:

1. State the mechanistic hypothesis, functional data source, disease-associated variant universe, and relevant contexts.
2. Establish the enhancer/silencer or other regulatory landscape with explicit element definitions and measured versus inferred activity.
3. Introduce the model or quantitative analysis, independent held-out performance, and baseline comparison.
4. Identify sequence determinants or regulatory grammar tied to measured function, then apply the model to candidate variants or elements.
5. Test context dependence, cross-species transfer, or prospective design when this is part of the selected route.
6. Present independent functional tests, exact-allele comparisons, and negative results before making broad mechanistic claims.
7. Work through at least two representative loci or elements, including close alternatives and target-gene evidence, and state the generalizable principle.

Change this order when the actual results support a clearer argument. Do not force a positive narrative. Null, conflicting, or underpowered findings are results and must be reported when they affect the conclusions.

When a prespecified hypothesis fails, the study still establishes something: identify those positive, evidence-bounded findings (for example colocalization-supported genes, nominated allele panels with comparators, validated model behaviour) and lead the title and the Abstract's final sentence with them, while reporting the failed hypothesis, its power and its bounds in the Abstract and Results. Do not hide the null result, and do not let it become the only message.

### Manuscript and Supplementary Materials files

Write the manuscript in one Markdown source built into a main manuscript file and a companion Supplementary Materials file (figs. S1 onward, then tables S1 onward), as in Science Advances. Include all text, publication-ready figures through embedded Markdown image links to final local assets, complete tables, legends, references, and disclosure statements in the following order. Generate and include every figure or table needed to support the scientific argument. Do not leave figure, table, caption, or citation placeholders in the completed manuscript.

#### 1. Title

Use a concise, informative, declarative title that names what the study discovered or nominates (for example the genes, loci, alleles or mechanism supported by the evidence) and the disease context. Prefer a positive, evidence-bounded statement over a title that only reports a failed test; use "nominate", "candidate" or "support" when evidence is predictive or associational. Italicize gene symbols. Avoid unsupported causal language, acronyms that are not widely recognized, and claims broader than the analyzed tissues, populations, or variant classes.

#### 2. Authors and affiliations

List authors, affiliations, corresponding-author information, and any consortium attribution exactly as provided by the investigators. Never infer authorship or author order. If these details have not been supplied, state that author and affiliation information remains to be provided.

#### 3. Abstract

Write one compact, unstructured paragraph containing:

- The disease problem and specific knowledge gap.
- The analytical or experimental approach.
- The principal quantitative results, including the most important counts and validation metric.
- The main biological interpretation.
- A final sentence stating the significance without exaggeration.

Use only results reported in the main manuscript. Do not cite references, introduce undefined specialist abbreviations, or make claims stronger than the evidence. Prefer exact values over adjectives such as "many," "strong," or "substantial." Keep it to at most 200 words, include the primary effect estimate and the power bound when the primary test is null, and end with the study's positive, bounded contribution.

#### 4. Introduction

Write a short, gap-driven Introduction that moves from disease importance to the specific unresolved regulatory-genomics problem. Closely follow the reference paper's progression:

1. Establish disease burden and defining biology.
2. Summarize genetic and genomic evidence relevant to the study.
3. Explain why noncoding variants and cell-type-specific regulation remain unresolved.
4. Describe the limitations of existing functional and computational approaches.
5. End with a precise statement of what this study did and what question it addresses.

Do not turn the Introduction into a comprehensive disease review. Include only background needed to understand the study. Reserve new results, detailed methods, and extended interpretation for their appropriate sections.

#### 5. Results

Use short, finding-based subheadings that state the subject and, when justified, the direction of the result. Each subsection should generally follow this sequence:

1. State the biological or analytical question.
2. Describe only the method needed to understand the result, referring detailed procedures to Materials and Methods.
3. Report the result with exact counts, denominators, effect sizes, uncertainty, and statistical support.
4. Refer to the relevant figure, table, or supplementary item at the point of the claim.
5. End with a restrained interpretation that advances the narrative.

Distinguish observations from interpretations. Use causal language only for results supported by a causal design or direct perturbation. Use "associated with," "consistent with," or "predicted to" for observational and model-derived findings. Report negative and sensitivity analyses where they affect confidence.

Make each Results sentence report a finding with its numbers; place procedural detail, software and parameters in Materials and Methods, and attach necessary qualifiers to the finding they bound rather than as standalone caveat sentences. Report every primary test with both model families when two exist, followed by its power or minimum detectable effect when it is null.

For every statistical comparison, provide the sample size, exact test, whether it was one- or two-sided, effect size, confidence interval when available, exact *P* value or a justified bound, and multiple-testing correction. Identify the unit of replication and distinguish biological from technical replicates. Report model discrimination and calibration on held-out data, including ROC AUC, PR AUC, class balance, decision thresholds, and uncertainty where available.

#### 6. Discussion

Open with the principal finding, not a repetition of the full Results. Then:

- Explain how the results address the original knowledge gap.
- Compare the findings with prior studies, including agreements and discrepancies.
- Interpret plausible mechanisms while labeling hypotheses as hypotheses.
- Discuss tissue, cell-type, ancestry, LD, annotation, model, and validation constraints.
- Explain the effect of missing information on the strength and generalizability of the conclusions.
- Identify the most informative next computational or experimental tests.
- End with a concise statement of what the study establishes.

Do not introduce unreported analyses in the Discussion. Do not claim clinical utility, causality, therapeutic relevance, or population generalizability unless directly supported.

#### 7. Materials and Methods

Provide enough information for an independent researcher to reproduce every reported result. Organize methods in the same order as the Results and include, as applicable:

- Study design, inclusion and exclusion criteria, and prespecified decisions.
- Literature-search databases, complete search dates, search strings, screening criteria, and evidence-selection rules.
- GWAS sources, cohort characteristics, ancestry definitions, phenotype definitions, sample sizes, genome builds, harmonization, quality control, significance thresholds, and LD reference panels.
- Locus construction, fine-mapping, variant annotation, coding and noncoding definitions, and variant-to-gene mapping.
- Epigenomic datasets, biosamples, accessions, cell types, peak-processing steps, coordinate conversion, blacklist filtering, and enhancer or silencer definitions.
- Model architecture, software and dependency versions, input encoding, train/validation/test partitions, random seeds, class balance, hyperparameters, checkpoints, thresholds, background variants, and evaluation metrics.
- TFBS, motif, target-gene, QTL, MPRA, conservation, repeat, population-genetic, and pathway analyses; distinguish exact-allele tests from interval overlap.
- Functional-assay design and analysis: construct or guide sequences, controls, randomization, donor and replicate structure, assay quality control, activity normalization, power or precision, tested nulls, primary endpoints, and validation independence.
- Statistical tests, assumptions, multiple-testing procedures, confidence intervals, missing-data handling, sensitivity analyses, and software versions.
- Data and code availability, including stable repository links and access restrictions.

State exact accession numbers, release versions, genome assemblies, coordinates, and dates. Cite both the original method or resource paper and the specific dataset release when appropriate. Do not use "default parameters" without naming the software version and all defaults that materially affect the result.

#### 8. Supplementary material

List the supplementary items under a `Supplementary Materials` heading in the main file ("This PDF file includes: figs. S1 to Sn; tables S1 to Sm") and place them in the companion Supplementary Materials file. Include supplementary methods, notes, figures, tables, robustness analyses, a result-to-claim map and extended results needed to evaluate reproducibility but not required for the main narrative. Number items sequentially as fig. S1 and table S1. Every supplementary item must be cited in the main or supplementary text and must have a self-contained legend.

#### 9. References

Use numbered citations in parentheses and in order of first appearance across the whole manuscript, including Materials and Methods, matching the reference paper's style, for example "(*1*, *2*)" or "(*3*–*5*)". Reuse the same number for repeated citations. Format entries in Science style: authors as initials then surname, title, italic abbreviated journal, volume, pages or article number, and year, for example "1. A. B. Author, C. D. Author, Title. *J. Abbr.* 12, 100–110 (2024)." The reference list must include complete author information as required by the target format, article title, journal or repository, volume, pages or article number, year, and DOI or stable URL when available.

Every citation must be real and verified against the primary source, PubMed, Crossref, the publisher record, or the official dataset or software repository. Open and inspect the source before citing it. Confirm that it supports the exact adjacent claim. Do not cite search-result snippets, AI-generated bibliographies, unverifiable references, or a review when the primary study is available and is the basis of the claim. Mark preprints clearly. For databases, software, and datasets, cite the authoritative resource and record the version or access date.

Perform a final bidirectional citation audit:

- Every factual claim requiring external support has an appropriate citation.
- Every in-text citation has exactly one matching reference-list entry.
- Every reference-list entry is cited in the text.
- Author names, title, year, journal, volume, pages or article number, and DOI agree with the authoritative record.
- Retractions, corrections, and expressions of concern have been checked.

#### 10. Acknowledgments and declarations

End the file with a Science Advances Acknowledgments block, using only information supplied or confirmed by the investigators: **Acknowledgments**, **Funding**, **Author contributions** (CRediT taxonomy), **Competing interests**, and **Data and materials availability** (primary data sources with accessions or URLs, access restrictions, derived data in the Supplementary Materials, and code repositories). Include ethics approval and informed consent when applicable.

List code repositories only after confirming that each is public: TREDNet v2 (https://github.com/tanoManzo/TREDNET_v2), the DeepFootprinting SHAP/TF-MoDISco/Tomtom pipeline (https://github.com/tanoManzo/FootPrinting) and the workflow framework (https://github.com/tanoManzo/disease-regulatory-genomics-workflow), plus the archived identifier of the deposited study scripts. Never redistribute restricted or controlled-access data (for example MetaBrain summary statistics or controlled-access GWAS data).

If required information is missing, say explicitly that it must be supplied before submission. Do not invent funding numbers, ethical approvals, contributions, conflicts, or repository links, and do not describe a planned deposit as completed.

### Figures and tables

Determine which findings require a visual display and generate the necessary publication-ready figures and tables directly from the verified outputs of Sections 0 through 6. Insert each item in the manuscript near its first citation. Use a figure when spatial, comparative, distributional, workflow, model-performance, or locus-level patterns are clearer visually. Use a table when readers need exact values, metadata, evidence provenance, or repeated-field comparisons. Describe simple findings in prose when a visual would add no information.

Figures and tables are required whenever they materially improve interpretation or make the reported evidence auditable. Their inclusion is determined by the results, not by a predetermined quota. If the available evidence does not justify a proposed item, omit it and state any consequential evidence gap in the text. Never fabricate values, extrapolate missing panels, or use a placeholder as though an analysis were complete.

Build the main figures around the manuscript's scientific argument rather than reproducing every analysis. The reference paper has eight main figures, and its major evidence roles suggest this adaptable plan:

- **Figure 1:** Define the variant or element universe, study design, measured data, model inputs, held-out partitions, and core performance. Put denominators and class balance beside performance metrics.
- **Figure 2:** Show regulatory-locus classes alongside matched expression, chromatin, or target-gene measurements in the relevant tissue. Make clear which links are measured and which are inferred.
- **Figure 3:** Show functional distinctions among element classes, including activity-linked sequence grammar, pathways, or context effects. Give matched backgrounds and corrected statistical support.
- **Figure 4:** Show cell-type, disease-state, or perturbation responses when these are measured. Separate compositional shifts from changes within a cell type.
- **Figure 5:** Show tested generalization across contexts or species, or prospective design outcomes when that is the chosen route. Omit this role when only one context is measured.
- **Figure 6:** Compare predictions with independent functional measurements at the exact tested element or allele, including effect sizes, directional agreement, uncertainty, and tested nulls.
- **Figures 7–8:** Work through at least two representative elements or loci from genetic association through regulatory activity, sequence mechanism, target gene, and close alternative variants. Label each evidence link and unresolved step.

Combine or omit roles when the evidence does not support them; do not create figures to meet a quota. A screen-centered manuscript may lead with assay design and functional calls instead of model performance. Tables should provide auditable numerical detail: dataset and assay metadata, GWAS/cohort summary when relevant, model splits and baselines, called elements with thresholds, independent validation coverage, and a prioritized candidate list with evidence levels. A table of proposed experiments is a future-work plan, not a Results table.

When the study used two model families, corroboration analyses and sensitivity analyses, give them main-figure roles rather than burying them: one figure for model performance of both families on held-out chromosomes and donors, prediction corroboration by independent assays and cross-family agreement; one figure for the primary tests repeated across families with locus-level intervals, direction against matched nulls and a power curve. Include a main table of the independent experimental datasets (Section 5) and a main table of candidate variants with comparators.

Design Fig. 1A as an illustrated study schematic in the reference paper's style, not a flowchart of text boxes: pastel rounded panels read left to right (inputs, models, allelic effects and prespecified tests, independent evidence, findings), each with a small drawn glyph (association plot, chromatin tracks, network diagram, REF/ALT sequence, effect lollipop, reporter construct, colocalization) and short labels with real numbers. Connect only adjacent panels with straight block arrows so no arrow crosses or overlaps a panel, glyph or label; feed parallel inputs into the next column with separate horizontal arrows. State in the legend that glyphs are illustrative and numbers come from the analyses. Keep schematic panels in a new figure file when earlier editions reuse the old one.

Generate figures from reproducible scripts and retain the plotted source data. Export final figures in a publication-suitable vector format when possible and in a high-resolution raster format when required. Embed the final local asset in Markdown using a relative path, for example `![Concise figure title](figures/figure_1.svg)`, followed immediately by its complete legend. Do not link to temporary files, interactive sessions, or absolute machine-specific paths.

Render compact tables directly in Markdown. For tables too large for the main text, include a concise main-text summary and place the complete table in the companion Supplementary Materials file. Preserve exact values and units; do not replace source values with manually rounded transcriptions unless the rounding rule is stated and applied consistently.

Every figure and table must be called out in numerical order and accompanied by a self-contained legend defining the population, tissue or cell type, data source, sample size, units, abbreviations, colors, error bars, statistical tests, and significance notation. Label axes and units, use colorblind-accessible palettes, show individual observations where feasible, and avoid plots that conceal distributions. Ensure that values remain legible at publication size. Cross-check every plotted or tabulated value against its source result before finalizing the manuscript.

#### Reference-based figure design and rendered-output review

Before generating figures, inventory the reference paper's actual main figures,
panel types, physical dimensions, typography, legends, and evidence progression.
Write a figure plan mapping each available scientific finding to its source data
and proposed display. Do not reduce the study to workflow diagrams and inventory
bars when verified outputs support informative allele-level, motif-level,
target-gene, tissue-specific, or representative-locus panels. Conversely, do not
match a figure quota by dividing one result into redundant figures or inventing
missing experiments. Document any evidence-dependent departures from the reference.

Design at the intended publication width, not at an oversized canvas that will
make labels unreadable after reduction. Record dimensions and minimum final-size
label fonts. Use distinct panel-label, title, axis, annotation, and legend areas;
long variant identifiers may use stable aliases only if the complete mapping is
provided. Do not place explanatory text over data or let legends obscure points,
intervals, or significance lines. Preserve meaningful uncertainty, denominators,
assessed alleles, coordinate systems, and the distinction between unavailable
coverage and a tested null result.

Inspect **every final figure and its actual rendered manuscript page** before
declaring the figure set complete. Review titles, labels, gene symbols, scientific
glyphs, clipping, collisions, legends, aspect ratios, rasterization, and caption
placement at the final size. A valid SVG or an overlap-free PNG does not establish
correct PDF rendering: verify that the renderer preserves fonts, sizes, weights,
and symbols rather than silently substituting defaults. Keep figures adjacent to
their complete legends and prefer vector embedding. Retain plotted source tables,
input provenance, reproducible scripts, and an explicit visual-review record.
Automated checks of asset presence, clipping, or typography are necessary where
feasible but do not replace this rendered-output review or scientific review.

Common rendering faults to check explicitly: crossing or overlapping arrows in
schematics; leader-line labels that collide with each other or with points; panel
letters overlapping axis labels; legends covering data; log axes without readable
tick labels; P values that underflow to zero (print a bound such as "< 10^-16");
Greek symbols missing from the PDF font (typeset them as math); unbreakable
accession lists that overflow table cells; and landscape tables that leave nearly
empty pages.

### Build, automated consistency checks and editions

Assemble the manuscript with a build script, not by hand. The build must (i) assert
that every reporting-stage number quoted in the prose matches its result table,
using exact formatted strings generated from the tables; (ii) assert counts of
figures, tables and citations and the absence of placeholders; (iii) number
citations and format references from the verified source table; (iv) typeset the
main file and Supplementary Materials and report layout warnings and page counts;
and (v) write a manifest with hashes of the source, builder, inputs, figures and
outputs. Build each new edition into new files and verify by hash that earlier
editions and any concurrently maintained versions are unchanged. Before editing a
shared manuscript directory, check file modification times for concurrent work by
another agent and do not overwrite it. Register each edition, its audit and any new
analysis in the study's results register and activity log.

### Writing style

Match the reference paper's concise, results-forward scientific style:

- Use direct topic sentences and place the main result early in each paragraph.
- Build paragraphs around one claim supported by quantitative evidence.
- Prefer active voice when the actor matters and precise passive voice when the procedure matters.
- Use consistent disease, locus, variant, model, tissue, and cell-type terminology.
- Define each abbreviation at first use and avoid unnecessary abbreviations.
- Report exact numbers and denominators before percentages.
- Separate data-supported conclusions from proposed mechanisms.
- Avoid promotional language, rhetorical questions, filler transitions, anthropomorphism, and repetition.
- Do not use em dashes. Follow the punctuation and typography rules stated earlier in this framework.
- Do not mention AI assistance or retain prompts, drafting notes, internal result identifiers, or workflow commentary in the submitted prose.

### Final scientific and editorial audit

Before declaring the manuscript complete, verify all of the following:

1. The title, abstract, Results, and Discussion express the same central claim at appropriate levels of detail.
2. Every quantitative statement matches its source output and uses the correct denominator.
3. Methods exist for every reported analysis, and every method has a corresponding result or clear purpose.
4. Claims distinguish association, prediction, functional evidence, and causality.
5. Population, ancestry, sex, age, tissue, cell-type, and disease-stage boundaries are reported and respected.
6. Figures, tables, legends, supplementary items, and references are complete and cited in order.
7. Multiple testing, uncertainty, class imbalance, data leakage, and sensitivity analyses are addressed where relevant.
8. Missing evidence is stated explicitly together with its impact; no values or citations are fabricated.
9. The prioritized variants can be traced to the evidence and thresholds generated in Sections 0 through 6.
10. The manuscript can be read as a self-contained scientific study without consulting the framework or analysis logs.
11. Language, nomenclature, gene symbols, variant identifiers, genome coordinates, abbreviations, and reference formatting are consistent.
12. The v5 readiness table links each C1-C10 component to actual evidence with the correct weights; the draft-only writing, claim-level novelty, journal gates, composite, and separate red-flag comments are audited. Proposed validation receives no completed-validation credit, and the C6 = 0 rule and its source inconsistency are disclosed when applicable.
13. Public code, shareable processed data, assay accession numbers, and restrictions on controlled data are stated accurately; private or restricted resources are not represented as public.
14. A final human review has confirmed authorship, interpretation, disclosures, ethics, data-sharing terms, and readiness for submission.
15. The manuscript follows the Science Advances layout of the reference paper (no Main or Extended Data sections), with an informative, positive and evidence-bounded title.
16. Every primary test is reported for each model family, with dependence-aware intervals and its power or minimum detectable effect when null.
17. Predicted elements are corroborated by assays independent of the training labels, and the independent experimental datasets are documented in a table with denominators, tested nulls, alignment, uncertainty and independence.
18. The build's automated number checks pass, the manifest is written, earlier editions are unchanged, and the rendered PDF pages, including Fig. 1A, have been inspected.
