#!/usr/bin/env python3
"""ALS-S2 step 02b: parse cases/controls and sex/age mentions from the
free-text sample-size fields of the studies table (the structured
cases/controls columns in the ancestries dump are empty for these studies).
"""
import pandas as pd, re

D = '/data/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/diseases/ALS/02_gwas'
s = pd.read_csv(f'{D}/results/ALS-S2-R001_studies_core.tsv', sep='\t', low_memory=False)

num = lambda x: int(x.replace(',', ''))
def parse_cc(txt):
    txt = str(txt)
    cases = sum(num(m) for m in re.findall(r'(\d[\d,]*)[^,;]*?\bcases\b', txt))
    ctrls = sum(num(m) for m in re.findall(r'(\d[\d,]*)[^,;]*?\bcontrols\b', txt))
    return cases, ctrls

s[['initial_cases','initial_controls']] = s['INITIAL SAMPLE SIZE'].apply(
    lambda t: pd.Series(parse_cc(t)))
s[['repl_cases','repl_controls']] = s['REPLICATION SAMPLE SIZE'].apply(
    lambda t: pd.Series(parse_cc(t)))

both = s['INITIAL SAMPLE SIZE'].fillna('') + ' ' + s['REPLICATION SAMPLE SIZE'].fillna('')
s['sex_mentioned'] = both.str.contains(r'\b(?:male|female|men|women)s?\b', case=False)
s['age_mentioned'] = both.str.contains(r'\bage', case=False)

dis = s[s.trait_subcategory == 'als_disease']
print('== ALS-S2-R001 addendum: cases/controls parsed from sample-size text ==')
print('(sums across studies; cohorts overlap heavily between consortium studies, DO NOT add up as unique individuals)')
print(f"als_disease studies: {len(dis)}")
print(f"initial: {dis.initial_cases.sum():,} cases / {dis.initial_controls.sum():,} controls")
print(f"replication: {dis.repl_cases.sum():,} cases / {dis.repl_controls.sum():,} controls")
print(f"largest single study by initial cases:")
top = dis.nlargest(3, 'initial_cases')[['STUDY ACCESSION','FIRST AUTHOR','DATE','initial_cases','initial_controls']]
print(top.to_string(index=False))
print(f"\nstudies with sex mentioned in cohort text: {s.sex_mentioned.sum()}/{len(s)}")
print(f"studies with age mentioned in cohort text: {s.age_mentioned.sum()}/{len(s)}")
s.to_csv(f'{D}/results/ALS-S2-R001_studies_core.tsv', sep='\t', index=False)
print('updated results/ALS-S2-R001_studies_core.tsv (added parsed case/control and sex/age flags)')
