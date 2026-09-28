#!/usr/bin/env python3
"""ALS-S3 step 01: overlap ALS candidate variants with known regulatory elements.
Variant set (GRCh38): 63 positioned genome-wide significant tags (ALS-S2 focal
set) plus 1000G LD proxies at r2>=0.8 (ALS-S2-R006, EUR union EAS).
Feature sets: ENCODE SCREEN cCREs V3, Ensembl Regulatory Build v116, FANTOM5
CAGE enhancers (hg38), UCSC RepeatMasker (hg38), ENCODE hg38 blacklist v2.
Tool: bedtools 2.31.1. Output: ALS-S3-R001 tables.
"""
import pandas as pd, subprocess, os

D2 = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
D3 = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/03_regulatory_landscape'
DATA = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/data'
BT = '/usr/local/apps/bedtools/2.31.1/bin/bedtools'

# --- variant BED ---
focal = pd.read_csv(f'{D2}/tmp/ld_focal.tsv', sep='\t')
prox = pd.read_csv(f'{D2}/results/ALS-S2-R006_ld_proxies.tsv', sep='\t')
rows = {}
for r in focal.itertuples():
    rows[r.snp] = (f'chr{r.chrom}', r.pos-1, r.pos, r.snp, 'tag')
for vid in prox.SNP_B.unique():
    p = str(vid).split(':')
    if len(p) != 4: continue
    c, pos, ref = p[0], int(p[1]), p[2]
    rows.setdefault(vid, (f'chr{c}', pos-1, pos-1+len(ref), vid, 'ld_proxy'))
v = pd.DataFrame(rows.values(), columns=['chrom','start','end','name','role']).sort_values(['chrom','start'])
bed = f'{D3}/tmp/candidates.bed'
v.to_csv(bed, sep='\t', header=False, index=False)
print(f'candidate variants: {len(v)} ({(v.role=="tag").sum()} tags, {(v.role=="ld_proxy").sum()} LD proxies)')

# --- feature files (normalized to chr-prefixed hg38 BED with a type column) ---
os.makedirs(f'{D3}/tmp', exist_ok=True)
feats = {}
feats['screen_ccre'] = (f'{DATA}/screen_ccres/GRCh38-cCREs.bed', None)  # col6 = group
# Ensembl GFF -> bed
ens_bed = f'{D3}/tmp/ensembl_reg.bed'
subprocess.run(f"zcat {DATA}/ensembl_regulation/Homo_sapiens.GRCh38.regulatory_features.v116.gff3.gz | "
               f"awk -F'\\t' '$1!~/^#/ {{print \"chr\"$1\"\\t\"$4-1\"\\t\"$5\"\\t.\\t.\\t\"$3}}' > {ens_bed}",
               shell=True, check=True)
feats['ensembl_reg'] = (ens_bed, None)
f5_bed = f'{D3}/tmp/fantom5.bed'
subprocess.run(f"zcat {DATA}/fantom5/F5.hg38.enhancers.bed.gz | cut -f1-3 | awk '{{print $0\"\\t.\\t.\\tenhancer\"}}' > {f5_bed}", shell=True, check=True)
feats['fantom5_enh'] = (f5_bed, None)
rmsk_bed = f'{D3}/tmp/rmsk.bed'
subprocess.run(f"zcat {DATA}/repeats/rmsk.txt.gz | awk -F'\\t' '{{print $6\"\\t\"$7\"\\t\"$8\"\\t.\\t.\\t\"$12}}' > {rmsk_bed}", shell=True, check=True)
feats['repeat'] = (rmsk_bed, None)
bl_bed = f'{D3}/tmp/blacklist.bed'
subprocess.run(f"zcat {DATA}/blacklist/hg38-blacklist.v2.bed.gz | awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t.\\t.\\tblacklist\"}}' > {bl_bed}", shell=True, check=True)
feats['blacklist'] = (bl_bed, None)

ann = v.set_index('name').copy()
for label, (path, _) in feats.items():
    out = subprocess.run([BT,'intersect','-a',bed,'-b',path,'-wa','-wb'],
                         capture_output=True, text=True, check=True).stdout
    hits = {}
    for line in out.strip().split('\n'):
        if not line: continue
        fl = line.split('\t')
        hits.setdefault(fl[3], set()).add(fl[-1])
    ann[label] = ann.index.map(lambda n: ';'.join(sorted(hits.get(n, []))))
    n_any = sum(1 for x in ann[label] if x)
    print(f'\n== {label}: {n_any}/{len(ann)} variants overlap ({n_any/len(ann):.1%}) ==')
    from collections import Counter
    c = Counter(t for x in ann[label] if x for t in x.split(';'))
    for k, n in c.most_common(12): print(f'  {k}: {n}')

ann.reset_index().to_csv(f'{D3}/results/ALS-S3-R001_variant_regulatory_overlap.tsv', sep='\t', index=False)
tags = ann[ann.role=='tag']
print('\n== tags only ==')
for label in feats:
    n = sum(1 for x in tags[label] if x)
    print(f'  {label}: {n}/{len(tags)} ({n/len(tags):.1%})')
print('\nsaved results/ALS-S3-R001_variant_regulatory_overlap.tsv')
