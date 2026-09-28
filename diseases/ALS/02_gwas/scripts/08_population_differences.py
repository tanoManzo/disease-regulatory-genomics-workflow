#!/usr/bin/env python3
"""ALS-S2 step 08: genome-wide significant gene sets by discovery-cohort
ancestry (European vs East Asian vs other/mixed), used in the population-
differences section of the Section 2 report. Log: logs/08_population_differences.log
"""
import pandas as pd, re
D = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
s = pd.read_csv(f'{D}/data/als_associations_single_variant.tsv', sep='\t', low_memory=False)
anc = pd.read_csv(f'{D}/data/als_ancestries.tsv', sep='\t', low_memory=False)
init = anc[anc.STAGE=='initial'].groupby('STUDY ACCESSION')['BROAD ANCESTRAL CATEGORY'].apply(
    lambda x: 'East Asian' if any('East Asian' in str(v) for v in x)
    else ('European' if all(str(v) in ('European','NR') for v in x) else 'other/mixed'))
gws = s[s.P <= 5e-8].copy()
gws['ancestry'] = gws['STUDY ACCESSION'].map(init)
def genes_of(g):
    g = str(g)
    return [] if g in ('nan','') else [x.strip() for x in re.split(r' - |, |; ', g) if x.strip()]
gx = gws.assign(genes=gws.MAPPED_GENE.map(genes_of)).explode('genes')
print('== GWS gene sets by discovery ancestry ==')
sets = {}
for a, grp in gx.groupby('ancestry'):
    sets[a] = set(grp.genes.dropna())
    print(f"\n{a}: {grp['STUDY ACCESSION'].nunique()} studies, {grp.SNPS.nunique()} GWS variants")
    print(' genes:', ', '.join(sorted(sets[a])))
if 'European' in sets and 'East Asian' in sets:
    print('\nshared genes (EUR & EAS):', ', '.join(sorted(sets['European'] & sets['East Asian'])) or 'none')
