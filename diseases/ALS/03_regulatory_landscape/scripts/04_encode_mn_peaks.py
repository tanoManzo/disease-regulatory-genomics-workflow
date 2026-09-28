#!/usr/bin/env python3
"""ALS-S3 step 04: download ENCODE processed peak files for ALS-donor iPSC
motor-neuron experiments (H3K27ac ChIP-seq and ATAC-seq), GRCh38.
Selection: released bed narrowPeak files, preferred_default where set,
output_type in (IDR thresholded peaks, replicated peaks, pseudoreplicated
peaks). Files land in 03_regulatory_landscape/data/encode_peaks/ with a
manifest TSV. ENCODE REST accessed via HPC proxy on run date.
Note: ENCODE has NO H3K27me3 for these biosamples (checked 2026-09-21), so
motor-neuron silencers cannot be defined from ENCODE data.
"""
import pandas as pd, requests, os, subprocess

D3 = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/03_regulatory_landscape'
OUT = f'{D3}/data/encode_peaks'; os.makedirs(OUT, exist_ok=True)
os.environ.setdefault('http_proxy','http://dtn26-e0:3128')
os.environ.setdefault('https_proxy','http://dtn26-e0:3128')

exps = pd.read_csv(f'{D3}/data/encode_als_motor_neuron_experiments.tsv', sep='\t')
sel = exps[((exps.assay=='Histone ChIP-seq') & (exps.target=='H3K27ac')) | (exps.assay=='ATAC-seq')]
print(f'selected experiments: {len(sel)} ({(sel.assay=="ATAC-seq").sum()} ATAC, {(sel.target=="H3K27ac").sum()} H3K27ac)')

GOOD_OUT = {'IDR thresholded peaks','replicated peaks','pseudoreplicated peaks','pseudo-replicated peaks'}
manifest = []
for r in sel.itertuples():
    j = requests.get(f'https://www.encodeproject.org/experiments/{r.accession}/?format=json', timeout=120).json()
    files = [f for f in j.get('files',[])
             if f.get('status')=='released' and f.get('assembly')=='GRCh38'
             and f.get('file_format')=='bed' and f.get('output_type') in GOOD_OUT]
    if not files:
        print(f'  {r.accession}: NO matching bed file'); continue
    pref = [f for f in files if f.get('preferred_default')]
    pick = (pref or sorted(files, key=lambda f: f.get('date_created',''), reverse=True))[0]
    url = 'https://www.encodeproject.org' + pick['href']
    dest = f"{OUT}/{r.accession}_{pick['accession']}_{(r.target if isinstance(r.target,str) and r.target else 'ATAC').replace(' ','')}.bed.gz"
    if not os.path.exists(dest):
        subprocess.run(['wget','-c','-q','--tries=3','--timeout=300','-O',dest,url], check=True)
    manifest.append({'experiment': r.accession, 'assay': r.assay, 'target': r.target if isinstance(r.target,str) and r.target else 'ATAC',
                     'file': pick['accession'], 'output_type': pick['output_type'],
                     'assembly': 'GRCh38', 'path': dest,
                     'biosample': r.biosample})
    print(f"  {r.accession} -> {pick['accession']} ({pick['output_type']})")
m = pd.DataFrame(manifest)
m.to_csv(f'{OUT}/manifest.tsv', sep='\t', index=False)
print(f'\ndownloaded {len(m)} peak files; manifest at data/encode_peaks/manifest.tsv')
