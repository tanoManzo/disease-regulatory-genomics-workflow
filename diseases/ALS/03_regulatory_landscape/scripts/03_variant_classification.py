#!/usr/bin/env python3
"""ALS-S3 step 03: classify ALS candidate variants (63 GWS tags + 1215 LD
proxies, GRCh38) as coding, enhancer-based, silencer-based, and/or repetitive.
coding    : overlaps a GENCODE v50 CDS exon.
enhancer  : overlaps H3K27ac peak (preliminary) / H3K27ac+DHS (refined) in the
            disease-relevant epigenome set (E069/E073 brain, E108 muscle).
silencer  : same logic with H3K27me3 (E069/E073/E107/E108).
repetitive: overlaps UCSC RepeatMasker element.
Percentages are reported per category (non-exclusive) and as a mutually
exclusive hierarchy coding > enhancer > silencer > repetitive-only > other.
"""
import pandas as pd, subprocess, os

D3 = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/03_regulatory_landscape'
DATA = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/data'
BT = '/usr/local/apps/bedtools/2.31.1/bin/bedtools'
TMP = f'{D3}/tmp'
bed = f'{TMP}/candidates.bed'

# CDS bed
cds = f'{TMP}/gencode_cds.bed'
if not os.path.exists(cds):
    subprocess.run(f"zcat {DATA}/gencode/gencode.v50.annotation.gtf.gz | "
                   f"awk -F'\\t' '$3==\"CDS\" {{print $1\"\\t\"$4-1\"\\t\"$5}}' | sort -k1,1 -k2,2n | {BT} merge -i - > {cds}",
                   shell=True, check=True)

def mk(label, cmd):
    out = f'{TMP}/class_{label}.bed'
    if not os.path.exists(out): subprocess.run(cmd + f' > {out}', shell=True, check=True)
    return out
enh_pre = mk('enh_pre', f'cat {TMP}/E069-H3K27ac.hg38.bed {TMP}/E073-H3K27ac.hg38.bed {TMP}/E108-H3K27ac.hg38.bed | cut -f1-3 | sort -k1,1 -k2,2n | {BT} merge -i -')
sil_pre = mk('sil_pre', f'cat {TMP}/E069-H3K27me3.hg38.bed {TMP}/E073-H3K27me3.hg38.bed {TMP}/E107-H3K27me3.hg38.bed {TMP}/E108-H3K27me3.hg38.bed | cut -f1-3 | sort -k1,1 -k2,2n | {BT} merge -i -')
enh_ref = mk('enh_ref', f'cat {TMP}/E069-H3K27ac.hg38.bed {TMP}/E073-H3K27ac.hg38.bed | {BT} intersect -u -a - -b {TMP}/DNase_brain.hg38.bed | cut -f1-3 > {TMP}/x1.bed; '
                        f'{BT} intersect -u -a {TMP}/E108-H3K27ac.hg38.bed -b {TMP}/DNase_muscle.hg38.bed | cut -f1-3 > {TMP}/x2.bed; '
                        f'cat {TMP}/x1.bed {TMP}/x2.bed | sort -k1,1 -k2,2n | {BT} merge -i -')
sil_ref = mk('sil_ref', f'cat {TMP}/E069-H3K27me3.hg38.bed {TMP}/E073-H3K27me3.hg38.bed | {BT} intersect -u -a - -b {TMP}/DNase_brain.hg38.bed | cut -f1-3 > {TMP}/y1.bed; '
                        f'cat {TMP}/E107-H3K27me3.hg38.bed {TMP}/E108-H3K27me3.hg38.bed | {BT} intersect -u -a - -b {TMP}/DNase_muscle.hg38.bed | cut -f1-3 > {TMP}/y2.bed; '
                        f'cat {TMP}/y1.bed {TMP}/y2.bed | sort -k1,1 -k2,2n | {BT} merge -i -')

v = pd.read_csv(bed, sep='\t', names=['chrom','start','end','name','role'])
def flag(path):
    out = subprocess.run([BT,'intersect','-u','-a',bed,'-b',path], capture_output=True, text=True, check=True).stdout
    hit = {l.split('\t')[3] for l in out.strip().split('\n') if l}
    return v.name.isin(hit)
v['coding'] = flag(cds)
v['enhancer_pre'] = flag(enh_pre); v['silencer_pre'] = flag(sil_pre)
v['enhancer_ref'] = flag(enh_ref); v['silencer_ref'] = flag(sil_ref)
v['repetitive'] = flag(f'{TMP}/rmsk.bed')

def report(sub, label):
    n = len(sub)
    print(f'\n== {label} (n={n}) ==')
    for c in ['coding','enhancer_pre','silencer_pre','enhancer_ref','silencer_ref','repetitive']:
        print(f'  {c}: {sub[c].sum()} ({sub[c].mean():.1%})')
    for suffix in ['pre','ref']:
        cls = pd.Series('other', index=sub.index)
        cls[sub['repetitive']] = 'repetitive_only'
        cls[sub[f'silencer_{suffix}']] = 'silencer'
        cls[sub[f'enhancer_{suffix}']] = 'enhancer'
        cls[sub['coding']] = 'coding'
        print(f'  exclusive hierarchy ({suffix}):', dict(cls.value_counts()))
report(v, 'all candidates (tags + LD proxies)')
report(v[v.role=='tag'], 'GWS tags only')
v.to_csv(f'{D3}/results/ALS-S3-R003_variant_classification.tsv', sep='\t', index=False)
print('\nsaved results/ALS-S3-R003_variant_classification.tsv')
