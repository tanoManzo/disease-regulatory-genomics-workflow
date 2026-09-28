#!/usr/bin/env python3
"""ALS-S2 step 04: top ALS GWAS genes.
Definition A (reproducibility): genes with genome-wide significant (P<=5e-8)
associations in >=2 distinct studies (STUDY ACCESSION) in the core ALS set.
Definition B (effect size): genome-wide significant associations ranked by
OR (>1 after flipping protective ORs) where OR reported.
Locus table: gene body (GENCODE v50 comprehensive, gene level) + left/right
flanking genes (nearest non-overlapping protein-coding or lncRNA gene) and
locus length (gene body span).
"""
import pandas as pd, re, os

D = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
GTF = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/data/gencode/gencode.v50.annotation.gtf.gz'

s = pd.read_csv(f'{D}/data/als_associations_single_variant.tsv', sep='\t', low_memory=False)
gws = s[s.P <= 5e-8].copy()

# explode mapped genes ("A - B" intergenic pairs, "A, B" lists)
def genes_of(g):
    g = str(g)
    if g in ('nan',''): return []
    return [x.strip() for x in re.split(r' - |, |; ', g) if x.strip() and x.strip()!='NA']
gws['genes'] = gws['MAPPED_GENE'].map(genes_of)
gx = gws.explode('genes').dropna(subset=['genes'])

rep = (gx.groupby('genes')
         .agg(n_studies=('STUDY ACCESSION','nunique'),
              n_assoc=('SNPS','size'),
              variants=('SNPS', lambda v: ','.join(sorted(set(v))[:6])),
              best_p=('P','min'))
         .sort_values(['n_studies','best_p'], ascending=[False,True]))
topA = rep[rep.n_studies >= 2]
print('== ALS-S2-R003A: genes genome-wide significant in >=2 studies ==')
print(topA.to_string())

# effect sizes: OR or BETA column; OR assumed when value > 0 and CI text lacks units
gws['effect'] = pd.to_numeric(gws['OR or BETA'], errors='coerce')
gws['ci_text'] = gws['95% CI (TEXT)'].astype(str)
gws['is_beta'] = gws.ci_text.str.contains(r'unit|increase|decrease|year|month|z score', case=False)
ors = gws[(~gws.is_beta) & gws.effect.notna()].copy()
ors['OR_ge1'] = ors.effect.where(ors.effect >= 1, 1/ors.effect)
topB = (ors.sort_values('OR_ge1', ascending=False)
        [['SNPS','MAPPED_GENE','effect','OR_ge1','P','FIRST AUTHOR','DATE','STUDY ACCESSION']]
        .drop_duplicates('SNPS').head(15))
print('\n== ALS-S2-R003B: largest ORs among genome-wide significant associations ==')
print(topB.to_string(index=False))

# --- GENCODE gene coordinates for locus table (top A genes + top B mapped genes) ---
want = set(topA.index) | set(g for gl in topB['MAPPED_GENE'].map(genes_of) for g in gl)
genes = []
import gzip
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line[0] == '#': continue
        fl = line.split('\t')
        if fl[2] != 'gene': continue
        attr = fl[8]
        name = re.search(r'gene_name "([^"]+)"', attr)
        gtype = re.search(r'gene_type "([^"]+)"', attr)
        genes.append((fl[0], int(fl[3]), int(fl[4]), name.group(1) if name else '',
                      gtype.group(1) if gtype else ''))
gdf = pd.DataFrame(genes, columns=['chrom','start','end','name','gene_type'])
gdf = gdf.drop_duplicates('name', keep='first')
main = gdf[gdf.gene_type.isin(['protein_coding','lncRNA'])].sort_values(['chrom','start']).reset_index(drop=True)

rows = []
for gname in sorted(want):
    hit = gdf[gdf.name == gname]
    if hit.empty:
        rows.append({'gene': gname, 'note': 'not in GENCODE v50 (may be alias)'}); continue
    h = hit.iloc[0]
    same = main[(main.chrom == h.chrom) & (main.name != gname)]
    left = same[same.end < h.start].nlargest(1, 'end')
    right = same[same.start > h.end].nsmallest(1, 'start')
    rows.append({'gene': gname, 'chrom': h.chrom, 'start': h.start, 'end': h.end,
                 'locus_length_bp': h.end - h.start + 1, 'gene_type': h.gene_type,
                 'left_flank': left.name.values[0] if len(left) else '',
                 'right_flank': right.name.values[0] if len(right) else ''})
loci = pd.DataFrame(rows)
print('\n== ALS-S2-R003C: locus table ==')
print(loci.to_string(index=False))

topA.reset_index().to_csv(f'{D}/results/ALS-S2-R003A_genes_multistudy.tsv', sep='\t', index=False)
topB.to_csv(f'{D}/results/ALS-S2-R003B_top_effect_sizes.tsv', sep='\t', index=False)
loci.to_csv(f'{D}/results/ALS-S2-R003C_locus_table.tsv', sep='\t', index=False)
print('\nsaved ALS-S2-R003{A,B,C} tables in results/')
