#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_PY="${SCRIPT_DIR}/hl_revision_pipeline.py"
WORKFLOW_SH="${SCRIPT_DIR}/run_hl_revision_workflow.sh"
TRAIN_PYTHON="/lustre/groups/crna01/projects/collabs/shared/cenvs/parnet_clean/bin/python"

usage() {
    cat <<'EOF'
Usage:
  run_all_hl_revision_workflow.sh [--train] [--umap] [--modisco-inputs] [--modisco-report --motif-db PATH] [--plots] [--overwrite] [--continue-on-error] [--cam-modes MODE[,MODE...]]

Examples:
  run_all_hl_revision_workflow.sh --train --umap --modisco-inputs
  run_all_hl_revision_workflow.sh --modisco-report --plots --motif-db /lustre/groups/crna01/projects/collabs/gabi/resubmission/data/external/pwms_all_motifs_ids.meme

Notes:
  - Runs over every entry in resubmission/configs/hl_revision_runs.json
  - Uses parnet_clean for config resolution and the existing single-run workflow wrapper for execution
  - With no stage flags, it regenerates all downstream outputs by default: UMAP, tf-modisco inputs, tf-modisco reports, and plot exports
  - By default the script stops on the first failed run; add --continue-on-error to keep going
EOF
}

DO_TRAIN=0
DO_UMAP=0
DO_MODISCO_INPUTS=0
DO_MODISCO_REPORT=0
DO_PLOTS=0
OVERWRITE=0
CONTINUE_ON_ERROR=0
MOTIF_DB="/lustre/groups/crna01/projects/collabs/gabi/resubmission/data/external/pwms_all_motifs_ids.meme"
CAM_MODES="final_logit_linearized"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --train)
            DO_TRAIN=1
            shift
            ;;
        --umap)
            DO_UMAP=1
            shift
            ;;
        --modisco-inputs)
            DO_MODISCO_INPUTS=1
            shift
            ;;
        --modisco-report)
            DO_MODISCO_REPORT=1
            shift
            ;;
        --motif-db)
            MOTIF_DB="$2"
            shift 2
            ;;
        --plots)
            DO_PLOTS=1
            shift
            ;;
        --overwrite)
            OVERWRITE=1
            shift
            ;;
        --continue-on-error)
            CONTINUE_ON_ERROR=1
            shift
            ;;
        --cam-modes)
            CAM_MODES="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ "${DO_TRAIN}" -eq 0 && "${DO_UMAP}" -eq 0 && "${DO_MODISCO_INPUTS}" -eq 0 && "${DO_MODISCO_REPORT}" -eq 0 && "${DO_PLOTS}" -eq 0 ]]; then
    DO_UMAP=1
    DO_MODISCO_INPUTS=1
    DO_MODISCO_REPORT=1
    DO_PLOTS=1
fi

if [[ "${DO_MODISCO_REPORT}" -eq 1 && -z "${MOTIF_DB}" ]]; then
    echo "--modisco-report requires --motif-db" >&2
    exit 1
fi

if [[ ! -x "${TRAIN_PYTHON}" ]]; then
    echo "Training python not found: ${TRAIN_PYTHON}" >&2
    exit 1
fi

mapfile -t RUN_NAMES < <("${TRAIN_PYTHON}" "${PIPELINE_PY}" list-runs | cut -f1)

if [[ "${#RUN_NAMES[@]}" -eq 0 ]]; then
    echo "No runs found in config." >&2
    exit 1
fi

FAILED_RUNS=()

for run_name in "${RUN_NAMES[@]}"; do
    [[ -z "${run_name}" ]] && continue
    echo "=== ${run_name} ==="

    cmd=("${WORKFLOW_SH}" --run "${run_name}")
    if [[ "${DO_TRAIN}" -eq 1 ]]; then
        cmd+=(--train)
    fi
    if [[ "${DO_UMAP}" -eq 1 ]]; then
        cmd+=(--umap)
    fi
    if [[ "${DO_MODISCO_INPUTS}" -eq 1 ]]; then
        cmd+=(--modisco-inputs)
    fi
    if [[ "${DO_MODISCO_REPORT}" -eq 1 ]]; then
        cmd+=(--modisco-report --motif-db "${MOTIF_DB}")
    fi
    if [[ "${DO_PLOTS}" -eq 1 ]]; then
        cmd+=(--plots)
    fi
    if [[ "${OVERWRITE}" -eq 1 ]]; then
        cmd+=(--overwrite)
    fi
    if [[ -n "${CAM_MODES}" ]]; then
        cmd+=(--cam-modes "${CAM_MODES}")
    fi

    if "${cmd[@]}"; then
        echo "Completed ${run_name}"
    else
        echo "Failed ${run_name}" >&2
        FAILED_RUNS+=("${run_name}")
        if [[ "${CONTINUE_ON_ERROR}" -ne 1 ]]; then
            exit 1
        fi
    fi
    echo

done

if [[ "${#FAILED_RUNS[@]}" -gt 0 ]]; then
    echo "Finished with failures: ${FAILED_RUNS[*]}" >&2
    exit 1
fi

echo "Finished all runs successfully."
