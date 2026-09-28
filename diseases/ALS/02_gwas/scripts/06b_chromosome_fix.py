#!/usr/bin/env python3
"""ALS-S2 step 06b: corrected chromosome distribution of GWS variants.
The catalog leaves CHR_ID empty for variants reported only as chr:pos, so the
chromosome is recovered from the variant name. Overwrites
results/ALS-S2-R005_chromosome_distribution.tsv (originally written by 06).
"""
import pandas as pd, re
D = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
s = pd.read_csv(f'{D}/data/als_associations_single_variant.tsv', sep='\t', low_memory=False)
gws = s[s.P <= 5e-8].drop_duplicates('SNPS').copy()
def chrom(r):
    if pd.notna(r['CHR_ID']): return str(r['CHR_ID']).replace('.0','')
    m = re.match(r'chr(\w+):', str(r['SNPS']))
    return m.group(1) if m else 'unknown'
gws['chrom_fixed'] = gws.apply(chrom, axis=1)
counts = gws['chrom_fixed'].value_counts()
print('== ALS-S2-R005B (corrected): GWS variant chromosome distribution ==')
print(counts.to_string())
counts.rename_axis('chrom').reset_index(name='n_gws_variants').to_csv(
    f'{D}/results/ALS-S2-R005_chromosome_distribution.tsv', sep='\t', index=False)
