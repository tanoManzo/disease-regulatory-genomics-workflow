#!/bin/bash
# One-time setup after cloning this repository on the NIH HPC (Biowulf).
#
# Links the large shared read-only assets (reference data, trained models,
# reference papers) into the clone so that all workflow paths resolve as if
# you were working in the original project folder. Nothing is copied and
# nothing is ever written into the shared folders: you only need read access
# to them. All work you generate (inputs you prepare, model outputs, results,
# logs) lives inside diseases/<your_disease>/ in this clone.
#
# Usage:   ./setup.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED=/vf/users/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex

for asset in data models reference_papers; do
    if [ ! -r "${SHARED}/${asset}" ]; then
        echo "ERROR: cannot read ${SHARED}/${asset}." >&2
        echo "Ask the project owner for read access to the shared folder." >&2
        exit 1
    fi
    ln -sfn "${SHARED}/${asset}" "${ROOT}/${asset}"
    echo "linked ${asset} -> ${SHARED}/${asset}"
done

# Recreate the writable working directories that git does not track.
for section in "${ROOT}"/diseases/*/0*_*/; do
    mkdir -p "${section}"/data "${section}"/results "${section}"/logs "${section}"/tmp
done
echo "created data/results/logs/tmp directories under diseases/"

echo
echo "Setup complete. To start a new disease study:"
echo "  1. Copy the section layout:  cp -r diseases/ALS diseases/<YOUR_DISEASE>  (then clear the ALS-specific TSV logs)"
echo "  2. Read frameworks/DiseaseFramework_SL091626.md and work through Sections 1-6."
echo "  3. Never write into data/ or models/ (they are read-only shared assets);"
echo "     all inputs and outputs belong in diseases/<YOUR_DISEASE>/0N_*/."
