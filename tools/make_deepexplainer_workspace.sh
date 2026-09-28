#!/bin/bash
# Create a writable DeepExplainer_TREDNet workspace inside your clone.
# submit2biowulf.sh resolves everything relative to its own location: it
# copies the code from backup/ into a timestamped output/ directory and
# submits the GPU job. So you need your own copy of the submit scripts
# (you edit their USER CONFIGURATION block anyway) next to a read-only
# symlink of backup/ (code + JASPAR MEME files), with a local output/.
#
# Usage:
#   ./tools/make_deepexplainer_workspace.sh diseases/<DISEASE>/04_modeling/deepexplainer
set -euo pipefail

WS=${1:?usage: make_deepexplainer_workspace.sh <workspace_dir>}
SHARED_DE=/vf/users/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/models/DeepExplainer_TREDNet

if [ ! -r "${SHARED_DE}/backup" ]; then
    echo "ERROR: cannot read ${SHARED_DE}; ask the owner for read access." >&2
    exit 1
fi

mkdir -p "${WS}/output"
ln -sfn "${SHARED_DE}/backup" "${WS}/backup"
cp -n "${SHARED_DE}/submit2biowulf.sh" "${SHARED_DE}/submit2biowulf_list.sh" "${WS}/"
chmod u+w "${WS}"/submit2biowulf*.sh

echo "DeepExplainer workspace ready: ${WS}"
echo
echo "Next steps:"
echo "  1. Edit ${WS}/submit2biowulf.sh USER CONFIGURATION:"
echo "     - EID / EXPERIMENT for your disease"
echo "     - TREDNET_PATH: keep the shared path for shared models, or for a"
echo "       model you trained yourself point TREDNET_MODEL_II_PATH and the"
echo "       INPUT_SEQS_* BEDs at your TREDNet workspace"
echo "       (diseases/<DISEASE>/04_modeling/trednet/{models_output,input_training_data})"
echo "     - the sbatch --qos flag, if you have a different GPU allocation"
echo "  2. bash ${WS}/submit2biowulf.sh"
echo "Results land in ${WS}/output/<EXPERIMENT>/<EID>_<timestamp>/ (SHAP plots,"
echo "TF-MoDISco motifs, motif.html). The shared conda env deepshap_env is"
echo "activated automatically by the job; read access is sufficient."
