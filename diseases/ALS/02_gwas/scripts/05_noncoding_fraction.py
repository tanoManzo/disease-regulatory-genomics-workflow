#!/usr/bin/env python3
"""ALS-S2 step 05: coding vs noncoding classification of core ALS GWAS variants.
Primary source: GWAS Catalog CONTEXT column (Ensembl most-severe consequence).
Coding classes: missense, synonymous, stop/start, frameshift, inframe indel,
coding_sequence, protein_altering. Splice donor/acceptor counted separately
(affect transcript, not protein-coding sequence directly).
Outputs ALS-S2-R004 tables: per-variant consequence, noncoding fraction overall
and at P<=5e-8, and per-gene coding/noncoding profile.
"""
import pandas as pd, re

D = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
s = pd.read_csv(f'{D}/data/als_associations_single_variant.tsv', sep='\t', low_memory=False)

CODING = {'missense_variant','synonymous_variant','stop_gained','stop_lost',
          'start_lost','frameshift_variant','inframe_insertion','inframe_deletion',
          'coding_sequence_variant','protein_altering_variant'}
SPLICE = {'splice_donor_variant','splice_acceptor_variant','splice_region_variant'}

v = (s.dropna(subset=['CONTEXT'])
      .drop_duplicates('SNPS')[['SNPS','CONTEXT','MAPPED_GENE','P']]
      .copy())
v['consequence'] = v.CONTEXT.astype(str).str.split(r'[;,]').str[0].str.strip()
def klass(c):
    if c in CODING: return 'coding'
    if c in SPLICE: return 'splice'
    return 'noncoding'
v['klass'] = v.consequence.map(klass)

print('== ALS-S2-R004: coding vs noncoding (unique variants with CONTEXT) ==')
print(f"variants with consequence annotation: {len(v)} of {s.SNPS.nunique()} unique variants")
print(v.klass.value_counts().to_string())
print(f"noncoding fraction: {(v.klass=='noncoding').mean():.1%}")
print('\nconsequence detail:'); print(v.consequence.value_counts().to_string())

gws = v[pd.to_numeric(v.P, errors='coerce') <= 5e-8]
print(f"\nat P<=5e-8 ({len(gws)} variants): noncoding fraction {(gws.klass=='noncoding').mean():.1%}")
print(gws.klass.value_counts().to_string())

# per-gene: fraction of GWAS genes with ONLY noncoding variants (GWS variants)
def genes_of(g):
    g = str(g)
    return [] if g in ('nan','') else [x.strip() for x in re.split(r' - |, |; ', g) if x.strip()]
gx = gws.assign(genes=gws.MAPPED_GENE.map(genes_of)).explode('genes').dropna(subset=['genes'])
prof = gx.groupby('genes')['klass'].agg(lambda k: 'only_noncoding' if set(k)=={'noncoding'}
                                        else ('has_coding' if 'coding' in set(k) else 'mixed_other'))
print(f"\n== genes (P<=5e-8): {len(prof)} total ==")
print(prof.value_counts().to_string())
print(f"fraction of GWS genes with only noncoding variants: {(prof=='only_noncoding').mean():.1%}")
print('\ngenes with coding GWS variants:', ', '.join(sorted(prof[prof=='has_coding'].index)))

v.to_csv(f'{D}/results/ALS-S2-R004_variant_consequences.tsv', sep='\t', index=False)
prof.reset_index().rename(columns={'klass':'profile'}).to_csv(
    f'{D}/results/ALS-S2-R004_gene_coding_profile.tsv', sep='\t', index=False)
print('\nsaved ALS-S2-R004 tables in results/')
