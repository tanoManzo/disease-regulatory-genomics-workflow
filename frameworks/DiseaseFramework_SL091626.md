# Framework of a foundational analysis prefacing a study of regulatory genomics of a particular disease.

## Introductory Query (if using an AI agent)
### Role and Final Goal
You are an autonomous research agent. Your final deliverable is a written scientific manuscript on the regulatory genomics of a specified disease. Every step you take exists to generate the evidence, tables, and figures needed to write that manuscript.

The disease under study is: [SPECIFY DISEASE]

### How to Work
Execute the attached framework in order (Sections 1 through 6). Create a subfolder per each section where you paste the results and the necessary data. Treat it as the analysis plan that produces the manuscript's results. Work through it methodically rather than jumping ahead.

Traceability. Save and log every step: the query run, the database or tool used, the parameters and thresholds chosen, the raw output, and the interpretation. Maintain a running record so any result in the final manuscript can be traced back to its source and reproduced. Assign each result a stable identifier and reference it when it appears in later sections or in the draft.

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
- `metabrain/` — complete MetaBrain GRCh38 cis-eQTL release for seven tissue/ancestry groups, with top effects, tabix indexes and the publication supplement containing trans-eQTL, interaction-eQTL and ALS enrichment tables. Reuse locally across diseases; redistribution is prohibited by the source terms.
- `fantom5/` — FANTOM5 CAGE enhancer atlas (hg38).
- `jaspar/` — JASPAR 2026 CORE vertebrate TF motifs (Sections 3-4 TFBS analysis).
- `blacklist/` — ENCODE hg38 blacklist v2 exclusion regions.
- `scripts/` — re-runnable download scripts (`wget -c` resumes); `logs/` — download logs.

Disease-specific inputs (e.g., full GWAS summary statistics for a given study, disease epigenomes) belong in `diseases/<disease>/0N_<section>/data/`, not here. Note for downloads on this HPC: compute nodes need the DTN proxy, `export http_proxy=http://dtnNN-e0:3128 https_proxy=$http_proxy` (NN = 17-27 or 30).


## Potentially Relevant Databases

- General: Psych-Encode, Cistrome, GEO, PRIDE, ClinVaR
- Cancer Specific: TCGA, Metabric, Cosmic

## 1. Background

### Epidemiology and heritability

- What are the incidence and prevalence of the disease? Stratify by age, sex, and geography (including ancestry/population where data exist).

- Is the disease heritable, and to what degree? Summarize evidence from twin and family studies and from SNP-based heritability estimates. Does the genetic risk architecture differ by ancestry?

### Genetics, pathways, and affected tissues

- What are the key genes associated with the disease?

- Which primary biological pathways are implicated?

- Which cell types and tissues are affected?

- Which ENCODE/NIH Roadmap cell lines are most relevant to those cell types and tissues?

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

- After expanding the tag SNPs into LD blocks, what are the counts of associated SNPs, indels, and translocations? Execute an LD analysis using 1000 genomes (or any relevant alternative) as reference.

- What are the top genes most strongly associated with this disease, by top genes we mean genes that have been reproduced in multiple GWAS studies at significant level? Now redefine top genes as those with largest effect sizes while being statistically significant. What are the effect sizes of these associations? Create a table listing these genes, their left- and right-flanking genes, and their locus lengths.

- What are the most extensively studied coding mutations associated with this disease?

- What fraction of GWAS genes have only non-coding mutations?

- What fraction of all GWAS variants is noncoding?

- Are GWAS gene loci (spanning the gene body + flanking regions) longer or shorter than the genome-wide average? Do they cluster on specific chromosomes? Are these genes evolutionarily well or poorly conserved? Have these genes been linked to any other diseases?

## 3. Regulatory Landscape

- Are there any genes whose expression is known to be altered in this disease, if so in which tissue/cell types?

- Are any of the GWAS variants located within known regulatory elements?

- List the studies profiling regulatory elements associated with this disease. Indicate which regulatory elements were identified in each study. What are the targets of these regulatory elements? Are any of these targets also genes whose expression is known to be altered in this disease? Are there any variants known to impact the activity of these regulatory elements?

- Using H3K27ac and H3K27me3 as preliminary marks of enhancers and silencers, respectively, count the number of enhancers and silencers within GWAS gene loci. For the previous query, redefine enhancers and silencers as regions wherein DHS (or ATAC-seq) peaks overlap with H3K27ac or H3k27me3 marks. Now count the number of enhancers and silencers using this updated definition.

- Focusing on GWAS SNPs/indels, calculate the percentage that are coding, enhancer-based, or silencer-based using this preliminary definition of regulatory elements. What percentage falls within repetitive elements (using UCSC definition)?

## 4. Modeling

- Train TREDNet (v2, two-phase transfer learning) enhancer and silencer models for the cell line(s)/cell type(s) associated with the disease. Inputs per model, built from cell-type-matched ENCODE data (hg38): positives = DNase or ATAC peaks overlapping H3K27ac (enhancer model) or H3K27me3 (silencer model), as 1 kb windows centered on the accessibility peak summit. Build controls from the reusable genome-wide pool of merged 1 kb DHS windows in `models/TREDNET_v2/input_training_data/control/allDHS.merge.nonPromoterExon`, which excludes promoters (GENCODE TSS +/- 2 kb) and exons. For each model, additionally remove control windows overlapping the ENCODE blacklist, the respective histone mark (H3K27ac or H3K27me3), or any positive window; restrict both classes to canonical chromosomes. Randomly sample controls with a recorded seed at an exact 1:2 positive-to-control ratio overall, and report positive/control counts for the training, validation, and test partitions. Train on chr1-6, chr10-22, chrX/Y; validate on chr7; test on chr8-9. Report test ROC AUC, PR AUC, and score thresholds at 10/5/3/1% FPR (used later as causal-variant thresholds). Consider modeling multiple cell lines/donors where appropriate and compare per-model performance. (Implementation: `models/TREDNET_v2/make_input_training_data.py` for inputs, using `--control-pool`, `--ratio 2`, and a recorded `--seed`; `TREDNet_v2.py` for training; `TREDNet_v2_inference.py` for scoring; e.g., trained ALS motor neuron enhancer model `MotorNeuron_Enhancer_DHS_x2`, ENCODE donor ENCDO689YDN ATAC ENCSR131HOY + H3K27ac ENCSR489LNU: AUC 0.782, AUPRC 0.652.)

- Using the TREDNet models trained in the previous step, re-count GWAS variants that occur within or affect predicted regulatory elements: a variant lies in a predicted enhancer/silencer if the model score of its surrounding region passes the FPR threshold chosen above, and its regulatory impact is quantified as the delta score between reference and alternate alleles (each scored in the variant-centered 2001 bp window). Establish a threshold for causal variants using these two quantities (e.g., region score above the 5% FPR cutoff and allele delta score in the top percentiles of matched background variants; state the exact definition used). Compare predicted causal regulatory variants with associated regulatory variants, coding variants, and all GWAS variants. Calculate the fraction of GWAS genes with at least one causal regulatory variant and the fraction with no regulatory variants. Classify gene loci into the following categories: coding + causal regulatory variants, coding only, causal regulatory only, and other. For the “other” loci, assess whether an alternative cell line might be necessary.

- Using DeepFootprinting (SHAP DeepExplainer applied to the trained TREDNet models), identify key transcription factors (TFs), histone modifications, and other genomic marks associated with enhancers and silencers: compute per-nucleotide SHAP contribution scores over the positive sequences, extract recurring high-importance subsequences as motifs with TF-MoDISco, and match them against JASPAR2024 with Tomtom to name the TFs. Calculate the fraction of predicted causal regulatory variants that fall within one of these key TF binding sites (TFBSs). What fraction of these variants disrupts the corresponding binding site? Do any create new TFBSs? (Implementation: `models/DeepExplainer_TREDNet/` — configure and submit with `submit2biowulf.sh` (single model) or `submit2biowulf_list.sh` (list of models); outputs SHAP scores, feature-importance summary plot, and a TF-MoDISco motif report `tfmodisco_output/motif.html` per run.)

- As a complementary approach, rank the key TFs by the number of predicted causal regulatory variants located within their binding sites. Determine whether any of the “other” GWAS variants overlap these binding sites outside of the predicted enhancers and silencers.

- Compare and contrast predicted causal regulatory variants in enhancers and silencers, as well as their associated TFBSs.

- How many of the predicted causal regulatory variants have been experimentally profiled? Are there experimental studies of noncoding variants that are not predicted to be causal? Why might our model be missing these variants, and how can they be recovered?

- Study the population genetics of candidate regulatory variants. How many are common, rare, or population-specific? How many are derived (DAF > 0.5)?

- Create THE LIST of rank-ordered, prioritized candidate regulatory variants for experimental follow-up.

- Map the identified enhancers/silencers to their target genes using proximity, re2g, HiChIP, etc.

## 5. Computational Validation

1. QTLs
   - Brain: MetaBrain, SingleBrain, PsychENCODE.
   - General: xQTL Serve, QTLbase, xQTL Atlas, GTEx.
   - Proteome QTLs: [ProteomeVariation](http://proteomevariation.org/).

2. MPRA datasets
   - [MPRAbase](http://mprabase.ucsf.edu/app/mprabase).

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

4. Other relevant datasets
   - BrainTF: ChIP-seq across cell types from different regions of the brain.
   - Genetic-disease catalogs: HGNC, DisGeNET, OMIM, Open Targets.
   - Variant–disease associations: [HuGE Calculator](https://hugeamp.org/hugecalculator.html), Open Targets.
   - Population genetics: [MAtCH](https://match.ctglab.nl/#/home), [Neale Lab UK Biobank LDSC heritability browser](https://nealelab.github.io/UKBB_ldsc/h2_browser.html).

### Shared local-data availability (updated 2026-09-25)

- Already local in `data/`: GTEx v10 cis-eQTLs, the complete MetaBrain 2021-07-23 cis-eQTL release and publication supplement, MPRAbase v4.9.3, all compact SingleBrain v2 top-association files, NIAGADS xQTL Atlas track metadata, consolidated BrainTF peak collections, HGNC, Open Targets, and Neale Lab UK Biobank LDSC topline results. The ALS-specific Project MinE 2021 GWAS, rare-burden, and SMR archive is stored in `data/als_gwas/project_mine_2021/`. Section 4 regulatory resources are also local.
- MetaBrain has restricted reuse terms: use the local copy for scientific research and education, do not redistribute it, and cite de Klein et al. (Nature Genetics, 2023; DOI 10.1038/s41588-023-01300-6). All 315 official manifest files passed MD5 verification.
- Not mirrored because access or licensing remains restricted: PsychENCODE primary data, OMIM, DisGeNET licensed releases, and participant-level biobank data. Obtain the required approval/account before use.
- Public query status for ALS: xQTL Atlas, QTLbase, ProteomeVariation, MAtCH, HuGE Calculator, and public PsychENCODE Phase I processed eQTLs are complete. xQTL Serve is unavailable and functionally superseded by xQTL Atlas.
- SingleBrain full association statistics (165.7 GB) are not mirrored; the complete top-association subset is local. Large NIAGADS/QTLbase bulk collections should be added only when the relevant QTL types and tissues have been selected.
- The authoritative path/access inventory is `data/validation_resources.tsv`; general SHA-256 checksums are in `data/validation_checksums.sha256`. Downloads are reproducible with `data/scripts/dl_validation_resources.sh` and `data/scripts/dl_metabrain.sh`; MetaBrain retains its official MD5 manifest and supplemental SHA-256 file inside its release directory.

### ALS execution status (2026-09-25)

Section 5 is complete for public data. Exact candidate queries, public processed QTLs, Project MinE GWAS-to-MetaBrain locus analysis, and applicable public portal audits are recorded in `diseases/ALS/05_computational_validation/`. Controlled/private biobanks, PsychENCODE primary data, and licensed OMIM/DisGeNET releases are excluded from this completion boundary.

## 6. Experimental and Biological Validation

- Prepare the list of the most promising candidate causal variants.

- Consider if TFs could be tested instead of the variants.

- Provide this data to collaborators capable of conducting MPRA, […] validation.

### ALS execution status (2026-09-25)

Section 6 design and handoff are complete in `diseases/ALS/06_experimental_validation/`: four ranked candidates, a TF-testing decision matrix, 16 GRCh38 MPRA inserts, and staged orthogonal-validation work packages. Wet-lab experiments have not been performed. External transfer or experimental initiation awaits collaborator selection.


## 7. Manuscript Writing and Assembly

### Objective

Write one complete, submission-ready scientific manuscript that converts the evidence generated in Sections 1 through 6 into a coherent study of disease-associated regulatory variation. The manuscript must not read like six independent reports or a chronological activity log. It must present a focused scientific argument in which each analysis motivates the next, the results support the conclusions, and the limitations are explicit.

Use `reference_papers/sciadv.adz3323.pdf` as the primary model for organization, pacing, quantitative detail, figure integration, and scientific tone. The style reference is:

> Huang D, Ovcharenko I. Silencer variants are key drivers of gene up-regulation in Alzheimer's disease. *Science Advances*. 2026;12:eadz3323. doi:10.1126/sciadv.adz3323.

Follow its observable writing conventions closely, but do not copy or lightly paraphrase its sentences. Do not import its Alzheimer's disease claims, results, citations, terminology, or interpretations into a manuscript about another disease. Cite the reference paper in the manuscript only if it is scientifically relevant, not merely because it was used as a writing model.

### Required inputs and readiness audit

Before drafting, inventory all outputs from Sections 1 through 6, including result files, tables, figures, statistical outputs, software logs, datasets, and verified literature references. Create an internal result-to-manuscript map containing:

- The stable result identifier.
- The exact quantitative finding and its units.
- The source file and analysis that produced it.
- The relevant sample size, genome build, tissue or cell type, ancestry, thresholds, and statistical test.
- The intended manuscript section, figure, or table.
- The source citation or dataset citation, where applicable.

Check that values reported in the text, figures, tables, legends, and supplementary material agree exactly. Resolve discrepancies from primary outputs rather than choosing the most convenient value.

Do not delay completion solely because an analysis or datum is unavailable. Produce a complete manuscript using the evidence that exists and identify every material gap explicitly. Use direct statements such as, "No reliable ancestry-stratified estimate was available," followed by the consequence for interpretation. Do not invent a value, conceal the gap, or leave unexplained `TODO`, `TBD`, or empty sections. If a missing analysis is essential to the central conclusion, state that the conclusion remains unestablished and specify the analysis needed to establish it.

### Central claim and narrative

Before writing prose, formulate one primary claim that is directly supported by the strongest results. Add two to four secondary claims that explain the genetic architecture, regulatory mechanisms, model performance, biological interpretation, or validation evidence. Every Results subsection must support one of these claims. Exclude technically correct analyses that do not advance the central narrative, or move them to the supplementary portion of the same file.

Organize the Results by scientific findings rather than by framework section number. A typical narrative is:

1. Define the disease-associated variant set, cohorts, loci, and relevant tissues or cell types.
2. Establish the noncoding and regulatory landscape of those loci.
3. Introduce and benchmark the disease-relevant enhancer and silencer models.
4. Identify and characterize candidate causal regulatory variants.
5. Connect variants to TF binding, target genes, pathways, cell types, and disease biology.
6. Present computational and, when available, experimental validation.
7. Prioritize candidates and explain their biological or translational implications.

Change this order when the actual results support a clearer argument. Do not force a positive narrative. Null, conflicting, or underpowered findings are results and must be reported when they affect the conclusions.

### Single-file manuscript structure

Write the complete manuscript in one Markdown file. Include all text, publication-ready figures through embedded Markdown image links to final local assets, complete tables, legends, references, supplementary material, and disclosure statements in the following order. Generate and include every figure or table needed to support the scientific argument. Do not leave figure, table, caption, or citation placeholders in the completed manuscript.

#### 1. Title

Use a concise, informative, and preferably declarative title that states the principal finding and disease context. Avoid unsupported causal language, acronyms that are not widely recognized, and claims broader than the analyzed tissues, populations, or variant classes.

#### 2. Authors and affiliations

List authors, affiliations, corresponding-author information, and any consortium attribution exactly as provided by the investigators. Never infer authorship or author order. If these details have not been supplied, state that author and affiliation information remains to be provided.

#### 3. Abstract

Write one compact, unstructured paragraph containing:

- The disease problem and specific knowledge gap.
- The analytical or experimental approach.
- The principal quantitative results, including the most important counts and validation metric.
- The main biological interpretation.
- A final sentence stating the significance without exaggeration.

Use only results reported in the main manuscript. Do not cite references, introduce undefined specialist abbreviations, or make claims stronger than the evidence. Prefer exact values over adjectives such as "many," "strong," or "substantial."

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
- TFBS, motif, target-gene, QTL, MPRA, conservation, repeat, population-genetic, and pathway analyses.
- Statistical tests, assumptions, multiple-testing procedures, confidence intervals, missing-data handling, sensitivity analyses, and software versions.
- Data and code availability, including stable repository links and access restrictions.

State exact accession numbers, release versions, genome assemblies, coordinates, and dates. Cite both the original method or resource paper and the specific dataset release when appropriate. Do not use "default parameters" without naming the software version and all defaults that materially affect the result.

#### 8. Supplementary material

Keep supplementary content in the same Markdown file under a clearly labeled `Supplementary Materials` heading. Include supplementary methods, notes, figures, tables, robustness analyses, and extended results needed to evaluate reproducibility but not required for the main narrative. Number items sequentially as Fig. S1, Table S1, and so forth. Every supplementary item must be cited in the main or supplementary text and must have a self-contained legend.

#### 9. References

Use numbered citations in order of first appearance, matching the reference paper's citation style. Reuse the same number for repeated citations. The reference list must include complete author information as required by the target format, article title, journal or repository, volume, pages or article number, year, and DOI or stable URL when available.

Every citation must be real and verified against the primary source, PubMed, Crossref, the publisher record, or the official dataset or software repository. Open and inspect the source before citing it. Confirm that it supports the exact adjacent claim. Do not cite search-result snippets, AI-generated bibliographies, unverifiable references, or a review when the primary study is available and is the basis of the claim. Mark preprints clearly. For databases, software, and datasets, cite the authoritative resource and record the version or access date.

Perform a final bidirectional citation audit:

- Every factual claim requiring external support has an appropriate citation.
- Every in-text citation has exactly one matching reference-list entry.
- Every reference-list entry is cited in the text.
- Author names, title, year, journal, volume, pages or article number, and DOI agree with the authoritative record.
- Retractions, corrections, and expressions of concern have been checked.

#### 10. Acknowledgments and declarations

End the file with the following labeled statements, using only information supplied or confirmed by the investigators:

- Acknowledgments.
- Funding.
- Author contributions using the CRediT taxonomy.
- Competing interests.
- Ethics approval and informed consent, when applicable.
- Data and materials availability.
- Code availability.

If required information is missing, say explicitly that it must be supplied before submission. Do not invent funding numbers, ethical approvals, contributions, conflicts, or repository links.

### Figures and tables

Determine which findings require a visual display and generate the necessary publication-ready figures and tables directly from the verified outputs of Sections 1 through 6. Insert each item in the manuscript near its first citation. Use a figure when spatial, comparative, distributional, workflow, model-performance, or locus-level patterns are clearer visually. Use a table when readers need exact values, metadata, evidence provenance, or repeated-field comparisons. Describe simple findings in prose when a visual would add no information.

Figures and tables are required whenever they materially improve interpretation or make the reported evidence auditable. Their inclusion is determined by the results, not by a predetermined quota. If the available evidence does not justify a proposed item, omit it and state any consequential evidence gap in the text. Never fabricate values, extrapolate missing panels, or use a placeholder as though an analysis were complete.

Build the main figures around the manuscript's scientific argument rather than reproducing every analysis. The exact number depends on the results, but the following pattern closely reflects the reference paper's evidence flow:

- **Figure 1:** Study design, cohort or GWAS input, variant filtering, locus definitions, and overall analysis workflow.
- **Figure 2:** Disease-relevant regulatory landscape, locus classes, cell-type or tissue context, and target-gene relationships.
- **Figure 3:** Enhancer and silencer model design, held-out performance, thresholds, and independent benchmarking.
- **Figure 4:** Candidate causal regulatory variants, allele-specific predictions, TFBS effects, and representative loci.
- **Figure 5:** Functional interpretation, gene expression, pathways, QTLs, chromatin contacts, and disease-state comparisons.
- **Figure 6:** Computational or experimental validation and the final prioritized candidate set.

Adapt, combine, or omit figures when the evidence does not support this layout. Do not create decorative panels or imply that a planned analysis was completed. Tables should provide auditable numerical detail, including a GWAS/cohort summary, top loci or genes, model datasets and performance, and the prioritized variant list with evidence scores.

Generate figures from reproducible scripts and retain the plotted source data. Export final figures in a publication-suitable vector format when possible and in a high-resolution raster format when required. Embed the final local asset in Markdown using a relative path, for example `![Concise figure title](figures/figure_1.svg)`, followed immediately by its complete legend. Do not link to temporary files, interactive sessions, or absolute machine-specific paths.

Render compact tables directly in Markdown. For tables too large for the main text, include a concise main-text summary and place the complete table in the supplementary portion of the same manuscript file. Preserve exact values and units; do not replace source values with manually rounded transcriptions unless the rounding rule is stated and applied consistently.

Every figure and table must be called out in numerical order and accompanied by a self-contained legend defining the population, tissue or cell type, data source, sample size, units, abbreviations, colors, error bars, statistical tests, and significance notation. Label axes and units, use colorblind-accessible palettes, show individual observations where feasible, and avoid plots that conceal distributions. Ensure that values remain legible at publication size. Cross-check every plotted or tabulated value against its source result before finalizing the manuscript.

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
9. The prioritized variants can be traced to the evidence and thresholds generated in Sections 1 through 6.
10. The manuscript can be read as a self-contained scientific study without consulting the framework or analysis logs.
11. Language, nomenclature, gene symbols, variant identifiers, genome coordinates, abbreviations, and reference formatting are consistent.
12. A final human review has confirmed authorship, interpretation, disclosures, ethics, data-sharing terms, and readiness for submission.
