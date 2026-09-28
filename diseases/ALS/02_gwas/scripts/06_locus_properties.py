#!/usr/bin/env python3
"""ALS-S2 step 06: properties of ALS GWAS gene loci.
- locus length vs genome-wide gene-length distribution (GENCODE v50,
  protein_coding + lncRNA gene bodies)
- chromosome distribution of genome-wide significant variants
- evolutionary conservation: mean phastCons100way / phyloP100way over gene
  bodies, ALS multi-study genes vs all protein-coding genes (UCSC
  bigWigAverageOverBed v2)
- other-disease links: other GWAS Catalog traits with P<=5e-8 in the same genes
"""
import pandas as pd, numpy as np, gzip, re, subprocess, os

D = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
DATA = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/data'
BWTOOL = '/usr/local/apps/ucsc/503/bin/x86_64/bigWigAverageOverBed'

loci = pd.read_csv(f'{D}/results/ALS-S2-R003C_locus_table.tsv', sep='\t')
topA = pd.read_csv(f'{D}/results/ALS-S2-R003A_genes_multistudy.tsv', sep='\t')
mgenes = set(topA['genes'])
loci_m = loci[loci.gene.isin(mgenes) & loci.chrom.notna()].copy()

# --- genome-wide gene lengths ---
genes = []
with gzip.open(f'{DATA}/gencode/gencode.v50.annotation.gtf.gz','rt') as f:
    for line in f:
        if line[0]=='#': continue
        fl = line.split('\t')
        if fl[2]!='gene': continue
        gt = re.search(r'gene_type "([^"]+)"', fl[8]).group(1)
        nm = re.search(r'gene_name "([^"]+)"', fl[8])
        genes.append((fl[0], int(fl[3]), int(fl[4]), nm.group(1) if nm else '', gt))
g = pd.DataFrame(genes, columns=['chrom','start','end','name','gene_type'])
g['length'] = g.end - g.start + 1
pc = g[g.gene_type=='protein_coding'].drop_duplicates('name')

print('== ALS-S2-R005A: locus length vs genome ==')
print(f"genome-wide protein-coding genes: n={len(pc)}, median length {pc.length.median():,.0f} bp, mean {pc.length.mean():,.0f} bp")
print(f"ALS multi-study gene loci (n={len(loci_m)}): median {loci_m.locus_length_bp.median():,.0f} bp, mean {loci_m.locus_length_bp.mean():,.0f} bp")
from scipy.stats import mannwhitneyu
pc_named = pc[~pc.name.isin(mgenes)]
u, p = mannwhitneyu(loci_m.locus_length_bp, pc_named.length, alternative='two-sided')
print(f"Mann-Whitney U vs all protein-coding: P = {p:.3g}")

# --- chromosome distribution of GWS variants ---
s = pd.read_csv(f'{D}/data/als_associations_single_variant.tsv', sep='\t', low_memory=False)
gws = s[s.P<=5e-8].drop_duplicates('SNPS')
print('\n== ALS-S2-R005B: chromosome distribution of GWS variants ==')
chrcounts = gws['CHR_ID'].astype(str).value_counts()
print(chrcounts.to_string())

# --- conservation ---
def bed_of(df, name_col, chrom_col='chrom', s_col='start', e_col='end'):
    return df.assign(score=0, strand='+')[[chrom_col, s_col, e_col, name_col]]
os.makedirs(f'{D}/tmp', exist_ok=True)
bed_als = f'{D}/tmp/als_genes.bed'; bed_all = f'{D}/tmp/all_pc_genes.bed'
loci_m[['chrom','start','end','gene']].assign(start=lambda d: d.start.astype(int)-1).to_csv(bed_als, sep='\t', header=False, index=False)
pc_named[['chrom','start','end','name']].assign(start=lambda d: d.start-1).drop_duplicates('name').to_csv(bed_all, sep='\t', header=False, index=False)
res = {}
for track in ['phastCons100way','phyloP100way']:
    bw = f'{DATA}/conservation/hg38.{track}.bw'
    for label, bed in [('als',bed_als), ('all',bed_all)]:
        out = f'{D}/tmp/{track}_{label}.tab'
        subprocess.run([BWTOOL, bw, bed, out], check=True)
        t = pd.read_csv(out, sep='\t', names=['name','size','covered','sum','mean0','mean'])
        res[(track,label)] = t
    a_m = res[(track,'als')]['mean']; g_m = res[(track,'all')]['mean']
    u, p = mannwhitneyu(a_m, g_m, alternative='two-sided')
    print(f"\n== ALS-S2-R005C: {track} gene-body mean ==")
    print(f"ALS genes: median {a_m.median():.4f} | all protein-coding: median {g_m.median():.4f} | MWU P = {p:.3g}")

cons = res[('phastCons100way','als')][['name','mean']].rename(columns={'mean':'phastCons_mean'}).merge(
       res[('phyloP100way','als')][['name','mean']].rename(columns={'mean':'phyloP_mean'}), on='name')
print(cons.sort_values('phastCons_mean', ascending=False).to_string(index=False))

# --- other diseases for ALS multi-study genes (GWAS Catalog cross-trait) ---
print('\n== ALS-S2-R005D: other GWS traits in ALS multi-study genes ==')
alla = pd.read_csv(f'{DATA}/gwas_catalog/gwas-catalog-download-associations-alt-full.tsv',
                   sep='\t', low_memory=False)
alla['P'] = pd.to_numeric(alla['P-VALUE'], errors='coerce')
allx = alla[alla.P<=5e-8].copy()
def genes_of(x):
    x = str(x)
    return [] if x in ('nan','') else [i.strip() for i in re.split(r' - |, |; ', x) if i.strip()]
allx['genes'] = allx.MAPPED_GENE.map(genes_of)
ax = allx.explode('genes')
rows = []
for gene in sorted(mgenes):
    hits = ax[(ax.genes==gene) & ~ax['MAPPED_TRAIT'].astype(str).str.contains('amyotrophic', case=False)]
    traits = hits['MAPPED_TRAIT'].value_counts().head(8)
    rows.append({'gene': gene, 'n_other_traits': hits['MAPPED_TRAIT'].nunique(),
                 'top_other_traits': '; '.join(traits.index[:6])})
ot = pd.DataFrame(rows)
print(ot.to_string(index=False))

loci_m.merge(cons, left_on='gene', right_on='name').to_csv(
    f'{D}/results/ALS-S2-R005_loci_conservation.tsv', sep='\t', index=False)
ot.to_csv(f'{D}/results/ALS-S2-R005_other_diseases.tsv', sep='\t', index=False)
chrcounts.rename_axis('chrom').reset_index(name='n_gws_variants').to_csv(
    f'{D}/results/ALS-S2-R005_chromosome_distribution.tsv', sep='\t', index=False)
print('\nsaved ALS-S2-R005 tables in results/')
