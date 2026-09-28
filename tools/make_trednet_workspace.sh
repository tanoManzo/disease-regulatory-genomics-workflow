#!/bin/bash
# Create a writable TREDNet v2 workspace inside your clone, wired to the
# shared read-only model assets. TREDNet_v2.py resolves model_phase_I/,
# input_training_data/ and models_output/ relative to the current directory,
# so training and inference must run from a directory YOU can write to.
# The shared folder stays untouched: your new input data, datasets,
# checkpoints and trained phase-II models all land inside the workspace.
#
# Usage:
#   ./tools/make_trednet_workspace.sh diseases/<DISEASE>/04_modeling/trednet
#
# Then, from the workspace (see the RUNBOOK.md it creates):
#   1. build inputs   -> input_training_data/
#   2. train          -> models_output/<EID>/
#   3. score variants -> your section's results/
set -euo pipefail

WS=${1:?usage: make_trednet_workspace.sh <workspace_dir>}
SHARED_TREDNET=/vf/users/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/models/TREDNET_v2

if [ ! -r "${SHARED_TREDNET}/model_phase_I" ]; then
    echo "ERROR: cannot read ${SHARED_TREDNET}; ask the owner for read access." >&2
    exit 1
fi

mkdir -p "${WS}"/input_training_data "${WS}"/models_output "${WS}"/logs
ln -sfn "${SHARED_TREDNET}/model_phase_I" "${WS}/model_phase_I"
ln -sfn "${SHARED_TREDNET}/fasta" "${WS}/fasta"

cat > "${WS}/RUNBOOK.md" <<EOF
# TREDNet workspace (created $(date +%F))

Run everything FROM THIS DIRECTORY. Code and the frozen phase-I model are
read from the shared folder; everything you generate is written here.

    SHARED=${SHARED_TREDNET}
    PY=\$SHARED/.venv/bin/python

## 1. Build training inputs for a new cell line / tissue
Peaks (DNase/ATAC + H3K27ac narrowPeak, hg38) go in your section's data/
directory first. Then:

    \$PY \$SHARED/make_input_training_data.py \\
      --eid <CellLine>_Enhancer_DHS_x2 \\
      --dnase <your_dnase.narrowPeak.gz> \\
      --h3k27ac <your_h3k27ac.narrowPeak.gz> \\
      --gencode-gtf ../../../../data/gencode/gencode.v50.annotation.gtf.gz \\
      --blacklist ../../../../data/blacklist/hg38-blacklist.v2.bed.gz \\
      --outdir input_training_data

## 2. Train phase II (GPU node; adapt \$SHARED/run_trednet_motorneuron.sh,
##    but cd HERE instead of the shared folder)

    export TREDNET_EID=<CellLine>_Enhancer_DHS_x2
    \$PY \$SHARED/TREDNet_v2.py        # writes models_output/\$TREDNET_EID/

## 3. Score sequences
Your new model:      --models-dir ./models_output
Existing shared models (e.g. MotorNeuron_Enhancer_DHS_x2):
                     --models-dir \$SHARED/models_output

    \$PY \$SHARED/TREDNet_v2_inference.py \\
      --eid <EID> --models-dir <see above> \\
      --seqs ../data/<variants.fa> --out ../results/<scores.tsv>
EOF

echo "TREDNet workspace ready: ${WS}"
echo "See ${WS}/RUNBOOK.md for the three-step recipe."
