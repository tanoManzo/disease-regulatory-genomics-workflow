#!/usr/bin/env python3
"""ALS-S2 step 07c: summarize LD expansion.
Input: tmp/ld_chr*_{EUR,EAS}.ld (plink r2>=0.8, +/-500kb windows) with variant
IDs CHROM:POS:REF:ALT. Classifies proxies as SNV (both alleles length 1) or
indel, counts per focal variant and per panel, and reports totals before and
after LD expansion. Output: results/ALS-S2-R006_ld_expansion.tsv (+ proxy list).
"""
import pandas as pd, glob, os

D = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
focal = pd.read_csv(f'{D}/tmp/ld_focal.tsv', sep='\t')

frames = []
for f in glob.glob(f'{D}/tmp/ld_chr*_*.ld'):
    pop = f.split('_')[-1].split('.')[0]
    t = pd.read_csv(f, sep=r'\s+')
    t['panel'] = pop
    frames.append(t)
ld = pd.concat(frames, ignore_index=True)

def vclass(vid):
    p = str(vid).split(':')
    return 'SNV' if len(p) == 4 and len(p[2]) == 1 and len(p[3]) == 1 else 'indel'
ld['proxy_class'] = ld['SNP_B'].map(vclass)
ld = ld[ld.SNP_A != ld.SNP_B]  # drop self-pairs

# map focal panel IDs back to catalog variant names by position
focal['panel_prefix'] = focal.chrom.astype(str) + ':' + focal.pos.astype(str) + ':'
id2snp = {}
for r in focal.itertuples():
    for vid in ld.SNP_A.unique():
        if str(vid).startswith(r.panel_prefix):
            id2snp[vid] = r.snp
ld['focal_snp'] = ld.SNP_A.map(id2snp)

print('== ALS-S2-R006: LD expansion (r2>=0.8, +/-500 kb, 1000G phase 3 GRCh38) ==')
n_focal_panel = ld.groupby('panel').SNP_A.nunique()
print('focal variants found in panel:'); print(n_focal_panel.to_string())
print(f"focal variants with position sought: {len(focal)} (of 71 GWS variants; 8 kgp IDs had no position)")

for pop in ['EUR','EAS']:
    sub = ld[ld.panel == pop]
    if sub.empty: continue
    uniq = sub.drop_duplicates('SNP_B')
    print(f"\n{pop}: {sub.SNP_A.nunique()} focal variants -> {sub.SNP_B.nunique()} unique proxies (incl. focal-as-proxy overlaps)")
    print(f"  proxy classes: {uniq.proxy_class.value_counts().to_dict()}")
    per = sub.groupby('focal_snp').SNP_B.nunique().sort_values(ascending=False)
    print(f"  median proxies per focal variant: {per.median():.0f}; top: {dict(per.head(5))}")

# combined unique proxy set (EUR union EAS)
uniq_all = ld.drop_duplicates('SNP_B')
print(f"\nunion over panels: {ld.SNP_B.nunique()} unique variants in LD blocks")
print(f"  classes: {uniq_all.proxy_class.value_counts().to_dict()}")
print("\nBefore expansion: 71 unique GWS variants (70 SNPs + 1 unresolved; 0 indels among classified).")
print("Translocations: none cataloged (GWAS Catalog does not curate translocations).")

ld.to_csv(f'{D}/results/ALS-S2-R006_ld_proxies.tsv', sep='\t', index=False)
summ = (ld.groupby(['panel','focal_snp'])
          .agg(n_proxies=('SNP_B','nunique'),
               n_snv=('proxy_class', lambda x: (x=='SNV').sum()),
               n_indel=('proxy_class', lambda x: (x=='indel').sum()),
               block_start=('BP_B','min'), block_end=('BP_B','max'))
          .reset_index())
summ['block_length_bp'] = summ.block_end - summ.block_start
summ.to_csv(f'{D}/results/ALS-S2-R006_ld_expansion.tsv', sep='\t', index=False)
print('\nsaved results/ALS-S2-R006_{ld_proxies,ld_expansion}.tsv')
