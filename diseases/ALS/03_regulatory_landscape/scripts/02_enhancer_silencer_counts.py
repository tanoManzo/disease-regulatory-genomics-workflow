#!/usr/bin/env python3
"""ALS-S3 step 02: enhancer and silencer counts within ALS GWAS gene loci.
Definitions per framework:
  preliminary: enhancer = H3K27ac narrowPeak; silencer = H3K27me3 narrowPeak.
  refined:     peak also overlaps a DNase (DHS) peak.
Epigenomes (Roadmap consolidated, hg19 lifted to hg38 with UCSC liftOver):
  marks: E069 (cingulate gyrus), E073 (DLPFC), E107/E108 (skeletal muscle M/F;
         E107 lacks H3K27ac).
  DNase pairing (no DNase exists for these four): brain marks use union of
  E081+E082 (fetal brain), muscle marks use union of E120+E121 (HSMM/myotube).
  Recorded as decision ALS-DEC-00x; fetal/cultured accessibility is a stated
  limitation.
GWAS gene loci: all genes with genome-wide significant variants (ALS-S2-R004
gene profile) with GENCODE v50 coordinates, gene body +/-100 kb flanks.
"""
import pandas as pd, subprocess, gzip, re, os

D2 = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
D3 = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/03_regulatory_landscape'
DATA = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/data'
BT = '/usr/local/apps/bedtools/2.31.1/bin/bedtools'
LIFT = '/usr/local/apps/ucsc/503/bin/x86_64/liftOver'
CHAIN = f'{DATA}/liftover/hg19ToHg38.over.chain.gz'
TMP = f'{D3}/tmp'; os.makedirs(TMP, exist_ok=True)

def lift(epi_mark):
    src = f'{DATA}/roadmap_peaks/{epi_mark}.narrowPeak.gz'
    if not os.path.exists(src):
        src = f'{DATA}/roadmap_peaks/{epi_mark}.macs2.narrowPeak.gz'
    out = f'{TMP}/{epi_mark}.hg38.bed'
    if os.path.exists(out): return out
    raw = f'{TMP}/{epi_mark}.hg19.bed'
    subprocess.run(f'zcat {src} | cut -f1-4 > {raw}', shell=True, check=True)
    subprocess.run([LIFT, raw, CHAIN, out, f'{out}.unmapped'], check=True)
    n_in = sum(1 for _ in open(raw)); n_out = sum(1 for _ in open(out))
    print(f'  lifted {epi_mark}: {n_out}/{n_in} peaks mapped')
    return out

marks = {'E069':'brain','E073':'brain','E107':'muscle','E108':'muscle'}
dnase_groups = {'brain':['E081-DNase','E082-DNase'], 'muscle':['E120-DNase','E121-DNase']}
print('lifting peaks hg19 -> hg38:')
files = {}
for e, grp in marks.items():
    for m in ['H3K27ac','H3K27me3']:
        if e=='E107' and m=='H3K27ac': continue  # not released
        files[(e,m)] = lift(f'{e}-{m}')
dnase = {}
for grp, es in dnase_groups.items():
    beds = [lift(x) for x in es]
    merged = f'{TMP}/DNase_{grp}.hg38.bed'
    subprocess.run(f'cat {beds[0]} {beds[1]} | sort -k1,1 -k2,2n | {BT} merge -i - > {merged}', shell=True, check=True)
    dnase[grp] = merged

# --- GWAS gene loci (gene body +/- 100kb) ---
prof = pd.read_csv(f'{D2}/results/ALS-S2-R004_gene_coding_profile.tsv', sep='\t')
want = set(prof['genes'])
genes = []
with gzip.open(f'{DATA}/gencode/gencode.v50.annotation.gtf.gz','rt') as f:
    for line in f:
        if line[0]=='#': continue
        fl = line.split('\t')
        if fl[2]!='gene': continue
        nm = re.search(r'gene_name "([^"]+)"', fl[8])
        if nm and nm.group(1) in want:
            genes.append((fl[0], max(0,int(fl[3])-100001), int(fl[4])+100000, nm.group(1)))
g = pd.DataFrame(genes, columns=['chrom','start','end','gene']).drop_duplicates('gene')
loci_bed = f'{TMP}/gws_gene_loci.bed'
g.sort_values(['chrom','start']).to_csv(loci_bed, sep='\t', header=False, index=False)
print(f'\nGWS gene loci: {len(g)} of {len(want)} genes have GENCODE v50 coordinates (+/-100kb flanks)')

# --- counts ---
def count_in_loci(peaks, refine_with=None):
    src = peaks
    if refine_with:
        src = f'{TMP}/tmp_refined.bed'
        subprocess.run(f'{BT} intersect -u -a {peaks} -b {refine_with} > {src}', shell=True, check=True)
    out = subprocess.run(f'{BT} intersect -u -a {src} -b {loci_bed} | wc -l',
                         shell=True, capture_output=True, text=True, check=True)
    total = sum(1 for _ in open(src))
    return int(out.stdout.strip()), total

rows = []
for (e, m), path in files.items():
    kind = 'enhancer(H3K27ac)' if m=='H3K27ac' else 'silencer(H3K27me3)'
    grp = marks[e]
    n_loci, n_tot = count_in_loci(path)
    rn_loci, rn_tot = count_in_loci(path, refine_with=dnase[grp])
    rows.append({'epigenome': e, 'tissue_group': grp, 'element': kind,
                 'peaks_total': n_tot, 'in_gws_loci_preliminary': n_loci,
                 'refined_total(DHS overlap)': rn_tot, 'in_gws_loci_refined': rn_loci})
res = pd.DataFrame(rows)
print('\n== ALS-S3-R002: enhancer/silencer counts in GWS gene loci ==')
print(res.to_string(index=False))
res.to_csv(f'{D3}/results/ALS-S3-R002_enhancer_silencer_counts.tsv', sep='\t', index=False)
g.to_csv(f'{D3}/results/ALS-S3-R002_gws_gene_loci.tsv', sep='\t', index=False)
print('\nsaved ALS-S3-R002 tables')
