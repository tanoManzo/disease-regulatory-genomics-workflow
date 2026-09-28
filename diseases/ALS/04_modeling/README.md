# Section 4: Modeling

Status: COMPLETE FOR THE MOTOR-NEURON ENHANCER ARM (2026-09-23).

The framework Section 4 workflow was run with TREDNet v2, REF/ALT scoring of 1,278 ALS candidates, SHAP DeepFootprinting on 1,000 positive sequences, TF-MoDISco/JASPAR2024 motif matching, EUR/EAS population annotation, and proximity + motor-neuron Hi-C target mapping.

Primary report: `results/section_4_modeling.md`  
Ranked experimental list: `results/ALS-S4-R009_THE_LIST.tsv`

Key result: held-out ROC AUC 0.780848 and PR AUC 0.651592; 30 predicted enhancer variants and four stringent causal regulatory candidates.

Hard limitation: no cell-matched motor-neuron H3K27me3 data exist in the ENCODE ALS collection, so a motor-neuron silencer model and enhancer-versus-silencer comparison are not estimable. This is recorded explicitly rather than substituted with a mismatched tissue.
