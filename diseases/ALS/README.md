# ALS regulatory genomics project

## Target disease

Amyotrophic lateral sclerosis (ALS). ALS is the primary phenotype. ALS-frontotemporal dementia (ALS-FTD) is included as contextual clinical and biological overlap but is not merged into ALS epidemiologic estimates.

## Framework

This workspace follows `frameworks/DiseaseFramework_SL091626.md` in order. Each numbered directory corresponds to one framework section and contains separate `data` and `results` directories. No scientific claim is considered usable in the manuscript until it has a stable result identifier and a verified source in the project logs.

## Analysis conventions

- Human genome build: GRCh38, with explicitly documented lift-over when a source uses GRCh37/hg19.
- Disease scope: sporadic and familial ALS, reported separately where evidence permits.
- Population scope: all ancestries, ages, sexes, and geographies, preserving source definitions and reporting missing strata.
- Evidence cutoff: the date on which each search is run, recorded in `activity_log.tsv`.
- Citations: numbered source identifiers in prose and DOI, PMID, accession, or stable URL in `sources.tsv`.

## Section 1 plan

1. Define ALS phenotype boundaries and terminology.
2. Collect incidence and prevalence estimates by geography, age, sex, and ancestry or population where reported.
3. Summarize familial, twin, and SNP-based heritability evidence, including ancestry limitations.
4. Compile replicated genes, pathways, vulnerable tissues and cell types, then map those cell types to relevant ENCODE and NIH Roadmap biosamples.
5. Review comorbidities, juvenile-onset disease, and clinical progression.
6. Document standard and emerging diagnostic approaches.
7. Build a treatment table covering status, molecular target, modality, mechanism, regulatory action, and pharmacogenomic evidence.
8. Synthesize open questions and evidence gaps without extrapolating beyond verified sources.

## Expected resources and access

Section 1 uses public resources and does not require private credentials. Expected sources include PubMed and primary literature, WHO and national registries, GeneReviews, ClinVar, Open Targets, ENCODE, NIH Roadmap Epigenomics, ClinicalTrials.gov, FDA, EMA, and other national regulatory agencies.

Later sections may require controlled-access authorization for All of Us, UK Biobank, or dbGaP. AlphaGenome or another selected sequence model may require model/API access and compute resources; those requirements will be logged before Section 4 begins.

## Current status

Sections 1-5 are complete for the agreed public-data scope. Section 6 is complete as an experimental-design and collaborator-handoff package as of 2026-09-25; wet-lab execution remains pending and no experimental result is claimed. Controlled/private biobanks, primary PsychENCODE data, and licensed catalogs remain excluded. Manuscript assembly has not started.
