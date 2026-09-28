#!/usr/bin/env python3
"""ALS-S4 prep: download chromatin loop calls (bedpe) for all ENCODE ALS-donor
motor-neuron in situ Hi-C experiments. Loops are used for enhancer/silencer
target-gene assignment in Sections 4-5. Full .hic contact matrices (2.5-8 GB
each) are deliberately NOT downloaded until the Section 4 modeling query is
finalized; accessions are recorded in the manifest for later retrieval.
"""
import pandas as pd, requests, os, subprocess

D3 = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/03_regulatory_landscape'
D4 = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/04_modeling'
OUT = f'{D4}/data/hic_loops'; os.makedirs(OUT, exist_ok=True)
os.environ.setdefault('http_proxy','http://dtn24-e0:3128')
os.environ.setdefault('https_proxy','http://dtn24-e0:3128')

exps = pd.read_csv(f'{D3}/data/encode_als_motor_neuron_experiments.tsv', sep='\t')
hic = exps[exps.assay == 'in situ Hi-C']
manifest = []
for acc in hic.accession:
    j = requests.get(f'https://www.encodeproject.org/experiments/{acc}/?format=json', timeout=120).json()
    for f in j.get('files', []):
        if f.get('status') != 'released': continue
        if f.get('output_type') == 'loops' and f.get('file_format') == 'bedpe':
            dest = f"{OUT}/{acc}_{f['accession']}_loops.bedpe.gz"
            if not os.path.exists(dest):
                subprocess.run(['wget','-c','-q','--tries=3','--timeout=120','-O',dest,
                                'https://www.encodeproject.org' + f['href']], check=True)
            manifest.append({'experiment': acc, 'file': f['accession'], 'kind': 'loops',
                             'preferred_default': bool(f.get('preferred_default')),
                             'size_MB': round((f.get('file_size') or 0)/1e6, 2), 'path': dest})
        elif f.get('output_type','').endswith('contact matrix'):
            manifest.append({'experiment': acc, 'file': f['accession'],
                             'kind': f['output_type'], 'preferred_default': bool(f.get('preferred_default')),
                             'size_MB': round((f.get('file_size') or 0)/1e6, 1), 'path': 'NOT_DOWNLOADED'})
m = pd.DataFrame(manifest)
m.to_csv(f'{OUT}/manifest.tsv', sep='\t', index=False)
dl = m[m.path != 'NOT_DOWNLOADED']
print(f'downloaded {len(dl)} loop files from {dl.experiment.nunique()} Hi-C experiments')
print(f'contact matrices recorded but deferred: {len(m)-len(dl)}')
print(m[m.preferred_default].to_string(index=False))
