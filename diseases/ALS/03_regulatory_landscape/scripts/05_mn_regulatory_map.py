#!/usr/bin/env python3
"""ALS-S3 step 05: motor-neuron regulatory map and candidate-variant overlap.
Enhancer definition in disease-relevant cells (ENCODE ALS-donor iPSC-derived
motor neurons, GRCh38): H3K27ac peak (union of 6 experiments) overlapping an
ATAC peak (union of 10 experiments). No H3K27me3 exists for these biosamples,
so a motor-neuron silencer map cannot be built from ENCODE (stated limit).
Outputs: MN enhancer BED, counts within GWS gene loci, and overlap of the
1,278 candidate variants (63 tags + LD proxies).
"""
import pandas as pd, subprocess, glob, os

D3 = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/03_regulatory_landscape'
BT = '/usr/local/apps/bedtools/2.31.1/bin/bedtools'
TMP = f'{D3}/tmp'; PK = f'{D3}/data/encode_peaks'

def union(pattern, out):
    files = sorted(glob.glob(f'{PK}/{pattern}'))
    subprocess.run(f'zcat {" ".join(files)} | cut -f1-3 | sort -k1,1 -k2,2n | {BT} merge -i - > {out}',
                   shell=True, check=True)
    return len(files), sum(1 for _ in open(out))

n_ac, n_ac_regions = union('*_H3K27ac.bed.gz', f'{TMP}/MN_H3K27ac_union.bed')
n_at, n_at_regions = union('*_ATAC.bed.gz', f'{TMP}/MN_ATAC_union.bed')
print(f'MN H3K27ac union: {n_ac_regions} regions from {n_ac} experiments')
print(f'MN ATAC union: {n_at_regions} regions from {n_at} experiments')

mn_enh = f'{TMP}/MN_enhancers.bed'
subprocess.run(f'{BT} intersect -u -a {TMP}/MN_H3K27ac_union.bed -b {TMP}/MN_ATAC_union.bed > {mn_enh}',
               shell=True, check=True)
n_enh = sum(1 for _ in open(mn_enh))
print(f'MN enhancers (H3K27ac ∩ ATAC): {n_enh} regions')

loci_bed = f'{TMP}/gws_gene_loci.bed'
for label, path in [('MN_H3K27ac_union', f'{TMP}/MN_H3K27ac_union.bed'),
                    ('MN_ATAC_union', f'{TMP}/MN_ATAC_union.bed'),
                    ('MN_enhancers', mn_enh)]:
    n = int(subprocess.run(f'{BT} intersect -u -a {path} -b {loci_bed} | wc -l',
                           shell=True, capture_output=True, text=True).stdout.strip())
    print(f'{label} within GWS gene loci (+/-100kb): {n}')

# variant overlap
bed = f'{TMP}/candidates.bed'
v = pd.read_csv(bed, sep='\t', names=['chrom','start','end','name','role'])
def flag(path):
    out = subprocess.run([BT,'intersect','-u','-a',bed,'-b',path], capture_output=True, text=True, check=True).stdout
    hit = {l.split('\t')[3] for l in out.strip().split('\n') if l}
    return v.name.isin(hit)
v['MN_H3K27ac'] = flag(f'{TMP}/MN_H3K27ac_union.bed')
v['MN_ATAC'] = flag(f'{TMP}/MN_ATAC_union.bed')
v['MN_enhancer'] = flag(mn_enh)
print('\n== candidate variants in motor-neuron regulatory regions ==')
for c in ['MN_H3K27ac','MN_ATAC','MN_enhancer']:
    print(f'  {c}: {v[c].sum()}/{len(v)} ({v[c].mean():.1%}) | tags: {v[v.role=="tag"][c].sum()}/63')
print('\nMN-enhancer variants by locus (top):')
hits = v[v.MN_enhancer]
print(hits[['chrom','start','name','role']].to_string(index=False))

subprocess.run(f'cp {mn_enh} {D3}/results/ALS-S3-R004_MN_enhancers.bed', shell=True, check=True)
v.to_csv(f'{D3}/results/ALS-S3-R004_variant_MN_overlap.tsv', sep='\t', index=False)
print('\nsaved ALS-S3-R004 outputs (MN enhancer BED + variant overlap)')
