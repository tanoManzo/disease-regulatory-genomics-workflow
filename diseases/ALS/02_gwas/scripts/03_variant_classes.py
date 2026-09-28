#!/usr/bin/env python3
"""ALS-S2 step 03: variant-level summary for core ALS associations.
- unique variants, variant class (SNP/indel/etc) via Ensembl REST (GRCh38),
- genome-wide significance filter, chromosome distribution, CNV flags.
Ensembl REST accessed 2026-09-21 via HPC proxy; classes cached to TSV so the
step is reproducible offline afterwards.
"""
import pandas as pd, requests, json, os, time

D = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
os.environ.setdefault('http_proxy','http://dtn20-e0:3128')
os.environ.setdefault('https_proxy','http://dtn20-e0:3128')

a = pd.read_csv(f'{D}/data/als_associations.tsv', sep='\t', low_memory=False)
a = a[a.als_category=='core_als'].copy()

# keep single-variant rows; multi-SNP haplotypes / interactions have ';' or ' x '
a['is_multi'] = a['SNPS'].astype(str).str.contains(r';| x ', regex=True)
print(f'association rows: {len(a)}, multi-variant rows: {a.is_multi.sum()}')
singles = a[~a.is_multi].copy()
variants = sorted(set(v for v in singles['SNPS'].astype(str) if v.startswith('rs')))
nonrs = sorted(set(singles['SNPS'].astype(str)) - set(variants))
print(f'unique single variants: {singles.SNPS.nunique()} ({len(variants)} rsIDs, {len(nonrs)} non-rs: {nonrs[:5]})')

cache_f = f'{D}/data/variant_classes_ensembl.tsv'
if os.path.exists(cache_f):
    vc = pd.read_csv(cache_f, sep='\t')
else:
    rows = []
    for i in range(0, len(variants), 200):
        chunk = variants[i:i+200]
        r = requests.post('https://rest.ensembl.org/variation/homo_sapiens',
                          headers={'Content-Type':'application/json'},
                          data=json.dumps({'ids': chunk}), timeout=120)
        r.raise_for_status()
        for rsid, d in r.json().items():
            m = d.get('mappings') or [{}]
            rows.append({'rsid': rsid, 'var_class': d.get('var_class'),
                         'minor_allele': d.get('minor_allele'), 'MAF': d.get('MAF'),
                         'chrom': m[0].get('seq_region_name'), 'pos_grch38': m[0].get('start'),
                         'alleles': m[0].get('allele_string')})
        print(f'  queried {i+len(chunk)}/{len(variants)}'); time.sleep(1)
    vc = pd.DataFrame(rows)
    missing = set(variants) - set(vc.rsid)
    for v in sorted(missing): rows.append({'rsid': v, 'var_class': 'NOT_FOUND'})
    vc = pd.DataFrame(rows)
    vc.to_csv(cache_f, sep='\t', index=False)

print('\n== ALS-S2-R002: variant classes (Ensembl 116, GRCh38) ==')
print(vc['var_class'].value_counts(dropna=False).to_string())

# genome-wide significant subset
singles['P'] = pd.to_numeric(singles['P-VALUE'], errors='coerce')
gws = singles[singles.P <= 5e-8]
gws_rs = set(gws['SNPS'])
print(f'\nassociations at P<=5e-8: {len(gws)} rows, {len(gws_rs)} unique variants')
vc_gws = vc[vc.rsid.isin(gws_rs)]
print('variant classes at P<=5e-8:')
print(vc_gws['var_class'].value_counts(dropna=False).to_string())
print('\nCNV-flagged association rows:', (a['CNV'].astype(str).str.upper()=='Y').sum())
print('translocations: the GWAS Catalog does not curate translocation associations; none present by design.')

singles.to_csv(f'{D}/data/als_associations_single_variant.tsv', sep='\t', index=False)
print('\nsaved: data/variant_classes_ensembl.tsv, data/als_associations_single_variant.tsv')
