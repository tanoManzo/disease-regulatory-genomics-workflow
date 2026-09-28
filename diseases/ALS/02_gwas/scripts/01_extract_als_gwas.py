#!/usr/bin/env python3
"""ALS-S2 step 01: extract ALS GWAS studies, associations, and ancestries
from the GWAS Catalog full dumps (release 2026-09-15, accessed 2026-09-21).

Trait definition (from EFO ontology, data/efo/efo.obo, accessed 2026-09-21):
MONDO:0004976 (amyotrophic lateral sclerosis) and all ontology descendants,
plus legacy EFO_0001357 (sporadic ALS) and EFO_0022918 (dominant ALS).
Categories:
  core_als     : ALS and pure-ALS subtypes (sporadic, familial, juvenile, typed)
  als_ftd      : ALS-FTD spectrum terms (kept separate per Section 1 scope)
  als_related  : associations whose reported trait mentions ALS but whose
                 mapped trait is a measurement/interaction (age of onset,
                 survival, pleiotropy partner diseases etc.) - not merged.
Run: /usr/local/Anaconda/envs/py3.11/bin/python3 01_extract_als_gwas.py
"""
import pandas as pd, re, os

DATA = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/data/gwas_catalog'
OUT = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas/data'
os.makedirs(OUT, exist_ok=True)

MONDO = 'http://purl.obolibrary.org/obo/MONDO_'
EFO = 'http://www.ebi.ac.uk/efo/EFO_'
core = {MONDO+m for m in ['0004976','0005144','0005145','0008780','0010459',
    '0011196','0011223','0011951','0012077','0012945','0013891','0014181',
    '0014531','0017593','0859529','0957538']} | {EFO+'0001357', EFO+'0022918'}
ftd = {MONDO+m for m in ['0007105','0013501','0014641','0030885']}
als_uris = core | ftd

def cat_of(uris):
    s = {u.strip() for u in str(uris).split(',')}
    if s & core: return 'core_als'
    if s & ftd: return 'als_ftd'
    return ''

# --- associations ---
a = pd.read_csv(f'{DATA}/gwas-catalog-download-associations-alt-full.tsv',
                sep='\t', low_memory=False)
a['als_category'] = a['MAPPED_TRAIT_URI'].map(cat_of)
mapped = a[a['als_category'] != ''].copy()

# reported-trait mentions of ALS not captured by the mapped URI (related traits)
rex = re.compile(r'amyotrophic', re.I)
rel = a[(a['als_category'] == '') & a['DISEASE/TRAIT'].astype(str).str.contains(rex)].copy()
rel['als_category'] = 'als_related'

mapped.to_csv(f'{OUT}/als_associations.tsv', sep='\t', index=False)
rel.to_csv(f'{OUT}/als_related_associations.tsv', sep='\t', index=False)

# --- studies ---
s = pd.read_csv(f'{DATA}/gwas-catalog-download-studies-v1.0.3.1.txt',
                sep='\t', low_memory=False)
s['als_category'] = s['MAPPED_TRAIT_URI'].map(cat_of)
smap = s[s['als_category'] != ''].copy()
srel = s[(s['als_category'] == '') & s['DISEASE/TRAIT'].astype(str).str.contains(rex)].copy()
srel['als_category'] = 'als_related'
studies = pd.concat([smap, srel])
studies.to_csv(f'{OUT}/als_studies.tsv', sep='\t', index=False)

# --- ancestries for those studies ---
anc = pd.read_csv(f'{DATA}/gwas-catalog-download-ancestries-v1.0.3.1.txt',
                  sep='\t', low_memory=False)
anc_als = anc[anc['STUDY ACCESSION'].isin(set(studies['STUDY ACCESSION']))].merge(
    studies[['STUDY ACCESSION','als_category']], on='STUDY ACCESSION', how='left')
anc_als.to_csv(f'{OUT}/als_ancestries.tsv', sep='\t', index=False)

# --- summary ---
print('== associations (mapped ALS traits) ==')
print(mapped['als_category'].value_counts().to_string())
print('unique variants (SNPS field), core_als:',
      mapped.loc[mapped.als_category=='core_als','SNPS'].nunique())
print('\n== related (reported-trait mentions, unmapped) ==')
print('rows:', len(rel), '| top reported traits:')
print(rel['DISEASE/TRAIT'].value_counts().head(10).to_string())
print('\n== studies ==')
print(studies['als_category'].value_counts().to_string())
print('unique publications (PubMed IDs):', studies['PUBMED ID'].nunique())
print('\n== mapped traits in core set ==')
print(mapped[mapped.als_category=='core_als']['MAPPED_TRAIT'].value_counts().to_string())
