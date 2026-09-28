#!/usr/bin/env python3
"""ALS-S2 step 07a: prepare LD expansion inputs.
Focal variants: unique core-ALS variants with P<=5e-8 and a GRCh38 position
(rsID positions from Ensembl 116 lookup; chr:pos IDs parsed directly;
kgp* array IDs have no public position and are excluded, counted).
LD reference: 1000 Genomes phase 3 GRCh38 (release 20190312), superpopulation
matched to discovery-cohort ancestry: EUR for all loci (discovery is
predominantly European), EAS additionally for loci discovered in East Asian
cohorts. Outputs: tmp/ld_focal.tsv, tmp/samples_EUR.txt, tmp/samples_EAS.txt,
tmp/regions_chr*.txt
"""
import pandas as pd, re, os

D = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
DATA = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/data'
os.makedirs(f'{D}/tmp', exist_ok=True)

s = pd.read_csv(f'{D}/data/als_associations_single_variant.tsv', sep='\t', low_memory=False)
vc = pd.read_csv(f'{D}/data/variant_classes_ensembl.tsv', sep='\t')
anc = pd.read_csv(f'{D}/data/als_ancestries.tsv', sep='\t', low_memory=False)

gws = s[s.P <= 5e-8].copy()
eas_studies = set(anc[(anc['STAGE']=='initial') &
                      (anc['BROAD ANCESTRAL CATEGORY'].astype(str).str.contains('East Asian'))]['STUDY ACCESSION'])
pos = {r.rsid: (str(r.chrom), int(r.pos_grch38)) for r in vc.itertuples()
       if pd.notna(r.pos_grch38) and str(r.chrom) in [str(c) for c in list(range(1,23))+['X']]}

rows, skipped = [], []
for snp, grp in gws.groupby('SNPS'):
    m = re.match(r'chr(\w+):(\d+)$', snp)
    if snp in pos: c, p = pos[snp]
    elif m: c, p = m.group(1), int(m.group(2))
    else: skipped.append(snp); continue
    pops = ['EUR'] + (['EAS'] if set(grp['STUDY ACCESSION']) & eas_studies else [])
    rows.append({'snp': snp, 'chrom': c, 'pos': p, 'pops': ','.join(pops),
                 'best_p': grp.P.min()})
f = pd.DataFrame(rows).sort_values(['chrom','pos'])
f.to_csv(f'{D}/tmp/ld_focal.tsv', sep='\t', index=False)
print(f'focal variants with position: {len(f)} | skipped (no position): {len(skipped)} -> {skipped}')
print('EAS-flagged focal variants:', (f.pops.str.contains('EAS')).sum())

panel = pd.read_csv(f'{DATA}/1000genomes/integrated_call_samples_v3.20130502.ALL.panel', sep='\t')
for sp in ['EUR','EAS']:
    ss = panel[panel.super_pop==sp]['sample']
    ss.to_csv(f'{D}/tmp/samples_{sp}.txt', index=False, header=False)
    print(f'{sp} samples: {len(ss)}')

for c, grp in f.groupby('chrom'):
    with open(f'{D}/tmp/regions_chr{c}.txt','w') as fh:
        for r in grp.itertuples():
            fh.write(f'chr{c}\t{max(0,r.pos-500000)}\t{r.pos+500000}\n')
print('wrote per-chromosome region files:', sorted(set(f.chrom)))
