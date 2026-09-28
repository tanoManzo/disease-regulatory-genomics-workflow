#!/usr/bin/env python3
"""ALS-S2 step 02: study-level summary - study counts, cohort ancestry
breakdown, sample sizes, cases/controls, sex/age availability.
Input: outputs of step 01. Output: results/ALS-S2-R001_* tables + stdout log.
"""
import pandas as pd, re, os

D = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
studies = pd.read_csv(f'{D}/data/als_studies.tsv', sep='\t', low_memory=False)
anc = pd.read_csv(f'{D}/data/als_ancestries.tsv', sep='\t', low_memory=False)
os.makedirs(f'{D}/results', exist_ok=True)

core_s = studies[studies.als_category == 'core_als'].copy()

# subcategorize compound core traits
def sub(t):
    t = str(t).lower()
    if 'age at onset' in t: return 'als_age_at_onset'
    if 'survival' in t: return 'als_survival'
    if 'drug use' in t: return 'als_drug_use_interaction'
    if 'multiple sclerosis' in t or 'frontotemporal' in t: return 'als_plus_other_disease'
    return 'als_disease'
core_s['trait_subcategory'] = core_s['MAPPED_TRAIT'].map(sub)

print('== ALS-S2-R001: GWAS study counts (GWAS Catalog release 2026-09-15) ==')
print('core ALS studies:', len(core_s), '| unique PubMed IDs:', core_s['PUBMED ID'].nunique())
print(core_s['trait_subcategory'].value_counts().to_string())
print('\nby year:'); print(core_s['DATE'].astype(str).str[:4].value_counts().sort_index().to_string())
print('\nfull summary stats available:', (core_s['FULL SUMMARY STATISTICS']=='yes').sum())

# --- ancestry breakdown (initial + replication stages, core ALS only) ---
anc_c = anc[anc.als_category == 'core_als'].copy()
anc_c['NUMBER OF INDIVIDUALS'] = pd.to_numeric(anc_c['NUMBER OF INDIVIDUALS'], errors='coerce')
for col in ['NUMBER OF CASES','NUMBER OF CONTROLS']:
    anc_c[col] = pd.to_numeric(anc_c[col], errors='coerce')

print('\n== ancestry breakdown, core ALS (individuals summed across studies; overlap between cohorts NOT removed) ==')
g = (anc_c.groupby(['STAGE','BROAD ANCESTRAL CATEGORY'])['NUMBER OF INDIVIDUALS']
     .agg(['count','sum']).rename(columns={'count':'n_cohort_entries','sum':'n_individuals'}))
print(g.to_string())

print('\n== cases/controls where reported (core ALS, all stages) ==')
cc = anc_c.groupby('BROAD ANCESTRAL CATEGORY')[['NUMBER OF CASES','NUMBER OF CONTROLS']].sum(min_count=1)
print(cc.to_string())

# sex/age: catalog has no structured sex/age columns; scan sample descriptions
desc = (anc_c['INITIAL SAMPLE DESCRIPTION'].fillna('') + ' ' +
        anc_c['REPLICATION SAMPLE DESCRIPTION'].fillna('') + ' ' +
        anc_c['SAMPLE DESCRIPTION'].fillna(''))
sex_mentions = desc.str.contains(r'\b(male|female|men|women)\b', case=False, regex=True)
age_mentions = desc.str.contains(r'\bage', case=False, regex=True)
print('\n== sex/age availability in cohort descriptions (core ALS cohort rows) ==')
print(f'rows: {len(anc_c)}, with sex mentions: {sex_mentions.sum()}, with age mentions: {age_mentions.sum()}')

# save tables
core_s.to_csv(f'{D}/results/ALS-S2-R001_studies_core.tsv', sep='\t', index=False)
g.reset_index().to_csv(f'{D}/results/ALS-S2-R001_ancestry_breakdown.tsv', sep='\t', index=False)
cc.reset_index().to_csv(f'{D}/results/ALS-S2-R001_cases_controls_by_ancestry.tsv', sep='\t', index=False)
print('\nsaved: results/ALS-S2-R001_{studies_core,ancestry_breakdown,cases_controls_by_ancestry}.tsv')
