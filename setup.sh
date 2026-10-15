#!/bin/bash
# One-time setup after cloning this repository on the NIH HPC (Biowulf).
#
# Links the large shared read-only assets (reference data and trained models)
# and the completed ALS worked example into the clone, so that every path in
# frameworks/DiseaseFramework.md resolves. Nothing is copied and nothing is
# ever written into the shared folders: you only need read access to them.
# All work you generate lives inside diseases/<YOUR_DISEASE>/ in this clone,
# which git ignores.
#
# Usage:   ./setup.sh                 link the shared assets
#          ./setup.sh <YOUR_DISEASE>  also create diseases/<YOUR_DISEASE>/
set -euo pipefail
shopt -s nullglob

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED=/vf/users/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex
EXAMPLE=ALS_ASTRA
SECTIONS=(00_study_design 01_background 02_gwas 03_regulatory_landscape 04_modeling
          05_computational_validation 06_experimental_validation 07_manuscript)

if [ "$(realpath "${ROOT}")" = "$(realpath "${SHARED}")" ]; then
    echo "This is the shared project folder itself; run setup.sh in your own clone." >&2
    exit 1
fi

for asset in data models; do
    if [ ! -r "${SHARED}/${asset}" ]; then
        echo "ERROR: cannot read ${SHARED}/${asset}." >&2
        echo "Ask the project owner for read access to the shared folder." >&2
        exit 1
    fi
    if [ -e "${ROOT}/${asset}" ] && [ ! -L "${ROOT}/${asset}" ]; then
        echo "ERROR: ${ROOT}/${asset} exists and is not a symlink; move it away first." >&2
        exit 1
    fi
    ln -sfn "${SHARED}/${asset}" "${ROOT}/${asset}"
    echo "linked ${asset} -> ${SHARED}/${asset}"
done

mkdir -p "${ROOT}/diseases"
if [ -r "${SHARED}/diseases/${EXAMPLE}" ] && [ ! -e "${ROOT}/diseases/${EXAMPLE}" ]; then
    ln -s "${SHARED}/diseases/${EXAMPLE}" "${ROOT}/diseases/${EXAMPLE}"
    echo "linked worked example diseases/${EXAMPLE} (read-only)"
fi

if [ $# -ge 1 ]; then
    NAME="$1"
    if [[ ! "${NAME}" =~ ^[A-Za-z0-9_-]+$ ]]; then
        echo "ERROR: the disease name may contain only letters, digits, '_' and '-'." >&2
        exit 1
    fi
    WS="${ROOT}/diseases/${NAME}"
    if [ -e "${WS}" ]; then
        echo "ERROR: ${WS} already exists." >&2
        exit 1
    fi
    for section in "${SECTIONS[@]}"; do
        mkdir -p "${WS}/${section}/scripts"
    done
    printf 'decision_id\tdate\tsection\tdecision\trationale\tstatus\tdecided_by\n' > "${WS}/decisions.tsv"
    printf 'action_id\tdate\tsection\taction\ttool_or_database\tparameters\toutput\tinterpretation\n' > "${WS}/activity_log.tsv"
    printf 'result_id\tdate\tsection\tdescription\tfile\tsource_ids\tnotes\n' > "${WS}/results_register.tsv"
    printf 'source_id\tkey\tdoi\tpmid\tauthors\ttitle\tjournal\tyear\tvolume\tpages\ttopic\tverified_via\tretraction_or_erratum\taccess_date\n' > "${WS}/sources.tsv"
    echo "created disease workspace diseases/${NAME}/"
fi

# Recreate the writable working directories of your own (non-linked) workspaces.
for disease in "${ROOT}"/diseases/*/; do
    [ -L "${disease%/}" ] && continue
    for section in "${disease}"0*_*/; do
        mkdir -p "${section}"data "${section}"results "${section}"logs "${section}"tmp
    done
done

echo
echo "Setup complete. Next steps:"
echo "  1. Create your workspace if you have not:  ./setup.sh <YOUR_DISEASE>"
echo "  2. Work through frameworks/DiseaseFramework.md, Sections 0 to 7, in order."
echo "  3. Never write into data/, models/ or diseases/${EXAMPLE}/ (read-only shared assets);"
echo "     all inputs and outputs belong in diseases/<YOUR_DISEASE>/0N_*/."
