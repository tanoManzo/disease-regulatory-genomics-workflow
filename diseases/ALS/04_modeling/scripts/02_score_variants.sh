#!/bin/bash
#SBATCH --job-name=als_trednet_score
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1,lscratch:20
#SBATCH --cpus-per-task=4
#SBATCH --mem=64g
#SBATCH --time=02:00:00
#SBATCH --output=../logs/02_score_variants_%j.log

set -euo pipefail

ROOT=/vf/users/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex
MODEL_DIR=${ROOT}/models/TREDNET_v2
SECTION=${ROOT}/diseases/ALS/04_modeling
EID=MotorNeuron_Enhancer_DHS_x2

cd "${MODEL_DIR}"
"${MODEL_DIR}/.venv/bin/python" TREDNet_v2_inference.py \
  --eid "${EID}" \
  --seqs "${SECTION}/data/ALS_variants_ref.fa" \
  --out "${SECTION}/results/ALS-S4-R002_ref_scores.tsv" \
  --chunk 256

"${MODEL_DIR}/.venv/bin/python" TREDNet_v2_inference.py \
  --eid "${EID}" \
  --seqs "${SECTION}/data/ALS_variants_alt.fa" \
  --out "${SECTION}/results/ALS-S4-R002_alt_scores.tsv" \
  --chunk 256
