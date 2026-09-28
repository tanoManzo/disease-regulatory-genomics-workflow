# Disease Regulatory Genomics Workflow

A reusable, agent-executable workflow for producing a scientific manuscript on
the regulatory genomics of a disease. The framework is disease-agnostic; the
completed ALS study is included as a worked example (scripts and READMEs; its
data and results stay on the shared filesystem).

## Layout

| Path | Contents | In git? |
|---|---|---|
| `frameworks/` | The master analysis framework (Sections 1-6) that drives the whole study | yes |
| `diseases/ALS/` | Worked example: per-section scripts, READMEs, and traceability logs | scripts + READMEs only |
| `data/` | ~127 GB of disease-agnostic reference data (GWAS Catalog, GENCODE, 1000G, GTEx, Roadmap, JASPAR, ...). See `data/README.md` for the inventory | no — symlink to shared folder |
| `models/` | Trained TREDNet v2 and DeepExplainer models (~2.3 TB) | no — symlink to shared folder |
| `reference_papers/` | Key publications | no — symlink to shared folder |

## Getting started (NIH Biowulf)

```bash
git clone <this-repo> my_disease_study
cd my_disease_study
./setup.sh
```

`setup.sh` symlinks `data/`, `models/`, and `reference_papers/` to the shared
read-only project folder
(`/vf/users/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/`) and
recreates the untracked `data/results/logs/tmp` directories in each section.
You need read access to that folder; nothing more. After that, every path in
the workflow resolves exactly as in the original project.

### Requirements

- **Dcode group membership.** `models/TREDNET_v2/` and
  `models/DeepExplainer_TREDNet/` are readable by group `Dcode` only;
  `setup.sh` checks this and tells you if access is missing.
- **Python for the analysis scripts** (`#!/usr/bin/env python3`): numpy,
  pandas, scipy, biopython, pysam. The shared env
  `/vf/users/Dcode/gaetano/conda/envs/generic_env` has all but pysam; either
  build your own small env or add pysam to a clone of it. The model
  pipelines need no setup: TREDNet uses its bundled `.venv` and DeepExplainer
  activates the shared `deepshap_env` automatically.
- **bcftools / plink** come from Biowulf modules (`module load bcftools plink`),
  as in the ALS scripts.

## The one rule that matters

**`data/` and `models/` are read-only shared assets. Never write there, and
never design a script that tries to.** All inputs you prepare and all outputs
you generate — including model inference inputs/outputs — belong in
`diseases/<YOUR_DISEASE>/0N_<section>/{data,results,logs,tmp}` inside your
clone. The ALS Section 4 scripts show the pattern: the model is loaded from
`models/TREDNET_v2/` but `--seqs` and `--out` point into the section's own
`data/` and `results/` directories.

## Training TREDNet for a new cell line or tissue

`TREDNet_v2.py` reads and writes relative to the current directory, so never
run it from inside the shared `models/TREDNET_v2/` folder (the first write
fails with Permission denied). Instead create a writable workspace in your
clone, wired to the shared phase-I model and genome FASTA:

```bash
./tools/make_trednet_workspace.sh diseases/<YOUR_DISEASE>/04_modeling/trednet
```

This creates `input_training_data/`, `models_output/`, and `logs/` locally,
symlinks the read-only `model_phase_I/` and `fasta/` from the shared folder,
and writes a `RUNBOOK.md` with the exact three commands: build inputs for
your cell line, train phase II, and score variants (with your new model or
any of the shared pre-trained ones, e.g. `MotorNeuron_Enhancer_DHS_x2`).
Your trained model stays in your clone under
`diseases/<YOUR_DISEASE>/04_modeling/trednet/models_output/<EID>/`.

## Running DeepExplainer (SHAP + TF-MoDISco) from your clone

`models/DeepExplainer_TREDNet/submit2biowulf.sh` writes its timestamped
output directories next to itself, so it must run from a copy in your clone,
not from the shared folder. Set one up with:

```bash
./tools/make_deepexplainer_workspace.sh diseases/<YOUR_DISEASE>/04_modeling/deepexplainer
```

This copies the two submit scripts (edit their USER CONFIGURATION block for
your EID, experiment name, and QOS), symlinks the read-only `backup/` (code
and JASPAR MEME databases), and creates a local `output/`. Point the model
and input BED paths at the shared TREDNet models or at the ones you trained
in your own TREDNet workspace. The job activates the shared `deepshap_env`
conda environment automatically.

## Starting a new disease study

1. Create `diseases/<YOUR_DISEASE>/` mirroring the ALS section layout
   (`01_background` ... `06_experimental_validation`, plus `manuscript/` and
   fresh `activity_log.tsv`, `decisions.tsv`, `results_register.tsv`,
   `sources.tsv`).
2. Give the framework (`frameworks/DiseaseFramework_SL091626.md`) to your AI
   agent with the disease name filled in, or work through it manually,
   Section 1 to Section 6 in order.
3. Follow the traceability conventions from the framework: stable result
   identifiers (`<DISEASE>-S<section>-R<nnn>`), logged queries and parameters,
   verified citations.

### Path convention for new scripts

The ALS scripts hardcode absolute paths to the original project folder (they
are the frozen record of a completed study; their references into shared
`data/` still work read-only from anywhere). For new scripts, resolve the
clone root from the script location so your study is portable:

```python
# diseases/<DISEASE>/0N_section/scripts/xx_step.py
from pathlib import Path
ROOT = Path(__file__).resolve().parents[4]   # clone root
DATA = ROOT / "data"                          # shared, read-only
SECTION = Path(__file__).resolve().parents[1] # this section (writable)
```

```bash
# diseases/<DISEASE>/0N_section/scripts/xx_step.sh
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
MODEL_DIR="${ROOT}/models/TREDNET_v2"          # shared, read-only
SECTION="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"  # writable
```

## HPC notes

- Compute nodes reach the internet through the DTN proxy:
  `export http_proxy=http://dtnNN-e0:3128 https_proxy=$http_proxy` (NN = 17-27 or 30).
- GPU jobs (TREDNet scoring) follow the SBATCH headers in
  `diseases/ALS/04_modeling/scripts/`.
- MetaBrain data in `data/metabrain/` may not be redistributed; use it in
  place only. Keep this repository private.
