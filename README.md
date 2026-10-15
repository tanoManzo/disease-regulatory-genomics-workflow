# Disease Regulatory Genomics Workflow

A reusable, agent-executable workflow for producing a scientific manuscript on
the regulatory genomics of a disease. The framework is disease-agnostic; the
completed ALS study is available read-only as a worked example (`setup.sh`
links it into your clone; it is not stored in git).

## Layout

| Path | Contents | In git? |
|---|---|---|
| `frameworks/` | The master analysis framework (`DiseaseFramework.md`, Sections 0-7) that drives the whole study | yes |
| `criteria/` | The manuscript evaluation rubric (`rubric_v5.md`) | yes |
| `reference_papers/` | The style reference paper (Huang and Ovcharenko, *Sci. Adv.* 2026, CC BY 4.0) | yes |
| `tools/` | Helpers that create writable TREDNet and DeepExplainer workspaces | yes |
| `diseases/<YOUR_DISEASE>/` | Your study: per-section scripts, data, results, logs and traceability TSVs | no — ignored; stays in your clone |
| `diseases/ALS_ASTRA/` | Worked example (scripts, results, manuscript editions) | no — read-only symlink |
| `data/` | ~127 GB of disease-agnostic reference data (GWAS Catalog, GENCODE, 1000G, GTEx, Roadmap, JASPAR, ...). See `data/README.md` for the inventory | no — symlink to shared folder |
| `models/` | TREDNet v2, DeepExplainer and AlphaGenome (~2.3 TB) | no — symlink to shared folder |

## Getting started (NIH Biowulf)

```bash
git clone https://github.com/tanoManzo/disease-regulatory-genomics-workflow.git my_disease_study
cd my_disease_study
./setup.sh <YOUR_DISEASE>
git switch <your-branch>
```

`setup.sh` symlinks `data/` and `models/` to the shared read-only project
folder (`/vf/users/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/`),
links the worked example as `diseases/ALS_ASTRA/`, and, when given a disease
name, creates `diseases/<YOUR_DISEASE>/` with the eight section folders
(`00_study_design` ... `07_manuscript`, each with `scripts/ data/ results/
logs/ tmp/`) and empty traceability logs (`decisions.tsv`, `activity_log.tsv`,
`results_register.tsv`, `sources.tsv`). Running it again is safe. You need
read access to the shared folder; nothing more.

### Branches

Each team member has a personal branch (`jaya`, `stas`, `ivan`, `elline`,
`alexandra`, `di`, `yash`); `dev` is the shared integration branch. Commit to
your branch, open a pull request into `dev` when work is ready to share, and
`dev` is merged into `main` after review. Pull `main` regularly.

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
  as in the worked-example scripts.
- **AlphaGenome** additionally needs a GPU allocation, your own Hugging Face
  account with the AlphaGenome model terms accepted, and a GRCh38 FASTA with
  its `.fai` index (see below).

## The one rule that matters

**`data/` and `models/` are read-only shared assets. Never write there, and
never design a script that tries to.** All inputs you prepare and all outputs
you generate — including model inference inputs/outputs — belong in
`diseases/<YOUR_DISEASE>/0N_<section>/{data,results,logs,tmp}` inside your
clone. The worked example's Section 4 scripts
(`diseases/ALS_ASTRA/04_modeling/scripts/`) show the pattern: the model is
loaded from `models/TREDNET_v2/` but inputs and outputs point into the
section's own `data/` and `results/` directories.

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

## Scoring variants with AlphaGenome

`models/AlphaGenome/` (see its `README.md` and framework Section 4) scores a
variant against 5,563 human tracks (RNA-seq, CAGE, DNase, ATAC, histone and TF
ChIP-seq, splicing) and reports a raw and a calibrated quantile score per gene
and track. Use the shared environment directly, on a GPU node, and write
outputs into your workspace:

```bash
AG=models/AlphaGenome
export ALPHAGENOME_FASTA=/path/to/GRCh38.fa
$AG/.venv/bin/hf auth login     # once; never put the token in files or scripts
$AG/.venv/bin/python $AG/score_variants.py --check-only --output-dir diseases/<YOUR_DISEASE>/04_modeling/results/ag_check
$AG/.venv/bin/python $AG/score_variants.py --chromosome <chr> --position <1-based> --ref <REF> --alt <ALT> \
    --output-dir diseases/<YOUR_DISEASE>/04_modeling/results/ag_<variant>
$AG/.venv/bin/python $AG/summarize_scores.py diseases/<YOUR_DISEASE>/04_modeling/results/ag_*/scores.tsv \
    --controls <controls>/*/scores.tsv --ontology <CURIEs> --output-dir <dir>
```

The runner scores one variant per call; for a variant set, write a batch
script in your workspace that loads the model once per GPU job.

## Starting a new disease study

1. Run `./setup.sh <YOUR_DISEASE>` to create `diseases/<YOUR_DISEASE>/` with
   the section layout and fresh traceability logs.
2. Give the framework (`frameworks/DiseaseFramework.md`) to your AI agent with
   the disease name filled in, or work through it manually, Section 0 to
   Section 7 in order. Worked examples in the framework point to
   `diseases/ALS_ASTRA/`; adapt them, do not reuse their disease-specific inputs.
3. Follow the traceability conventions from the framework: stable result
   identifiers (`<DISEASE>-S<section>-R<nnn>`), logged queries and parameters,
   verified citations.

### Path convention for new scripts

Some worked-example scripts hardcode absolute paths to the original project
folder (they are the frozen record of a completed study; their references into
shared `data/` still work read-only from anywhere). For new scripts, resolve
the clone root from the script location so your study is portable:

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
- GPU jobs (TREDNet training and scoring) follow the SBATCH headers in
  `diseases/ALS_ASTRA/04_modeling/scripts/*.sbatch`.
- MetaBrain data in `data/metabrain/` may not be redistributed; use it in
  place only.
- **This repository is public.** Never commit credentials (including Hugging
  Face tokens), restricted or controlled-access data, unpublished manuscripts,
  or anything under `diseases/`.
