#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PIPELINE_PY="${SCRIPT_DIR}/hl_revision_pipeline.py"
PLOT_EXPORT_PY="${SCRIPT_DIR}/export_hl_revision_plots.py"
TRAIN_PYTHON="/lustre/groups/crna01/projects/collabs/shared/cenvs/parnet_clean/bin/python"

usage() {
    cat <<'EOF'
Usage:
  run_hl_revision_workflow.sh --run RUN [--train] [--umap] [--modisco-inputs] [--modisco-report --motif-db PATH] [--plots] [--fold N] [--checkpoint PATH] [--cam-modes MODE[,MODE...]]

Examples:
  run_hl_revision_workflow.sh --run hl_revised_20percent --umap
  run_hl_revision_workflow.sh --run hl_revised_raw_gabistab --modisco-inputs --modisco-report --plots --motif-db /path/to/pwms_all_motifs_ids.meme

Environment split:
  - parnet_clean: validation, training, UMAP, CAM/tf-modisco input export
  - tfmodisco: motif discovery and report only
EOF
}

RUN_NAME=""
DO_TRAIN=0
DO_UMAP=0
DO_MODISCO_INPUTS=0
DO_MODISCO_REPORT=0
DO_PLOTS=0
FOLD=""
CHECKPOINT=""
MOTIF_DB=""
OVERWRITE=0
CAM_MODES="final_logit_linearized"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run)
            RUN_NAME="$2"
            shift 2
            ;;
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
        --plots)
            DO_PLOTS=1
            shift
            ;;
        --motif-db)
            MOTIF_DB="$2"
            shift 2
            ;;
        --fold)
            FOLD="$2"
            shift 2
            ;;
        --checkpoint)
            CHECKPOINT="$2"
            shift 2
            ;;
        --overwrite)
            OVERWRITE=1
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

if [[ -z "${RUN_NAME}" ]]; then
    usage >&2
    exit 1
fi

if [[ ! -x "${TRAIN_PYTHON}" ]]; then
    echo "Training python not found: ${TRAIN_PYTHON}" >&2
    exit 1
fi

"${TRAIN_PYTHON}" "${PIPELINE_PY}" validate --run "${RUN_NAME}"

COMMON_ARGS=(--run "${RUN_NAME}")
if [[ -n "${FOLD}" ]]; then
    COMMON_ARGS+=(--fold "${FOLD}")
fi
if [[ -n "${CHECKPOINT}" ]]; then
    COMMON_ARGS+=(--checkpoint "${CHECKPOINT}")
fi

if [[ "${DO_TRAIN}" -eq 1 ]]; then
    TRAIN_ARGS=("${COMMON_ARGS[@]}")
    if [[ -n "${FOLD}" ]]; then
        TRAIN_ARGS=(--run "${RUN_NAME}" --folds "${FOLD}")
    fi
    if [[ "${OVERWRITE}" -eq 1 ]]; then
        TRAIN_ARGS+=(--overwrite)
    fi
    "${TRAIN_PYTHON}" "${PIPELINE_PY}" train "${TRAIN_ARGS[@]}"
fi

if [[ "${DO_UMAP}" -eq 1 ]]; then
    "${TRAIN_PYTHON}" "${PIPELINE_PY}" export-umap "${COMMON_ARGS[@]}"
fi

if [[ "${DO_MODISCO_INPUTS}" -eq 1 ]]; then
    "${TRAIN_PYTHON}" "${PIPELINE_PY}" export-modisco-inputs "${COMMON_ARGS[@]}" --cam-modes "${CAM_MODES}"
fi

if [[ "${DO_MODISCO_REPORT}" -eq 1 ]]; then
    IFS=',' read -r -a CAM_MODE_LIST <<< "${CAM_MODES}"
    for cam_mode in "${CAM_MODE_LIST[@]}"; do
        cam_mode="$(echo "${cam_mode}" | xargs)"
        [[ -z "${cam_mode}" ]] && continue
        MODISCO_ARGS=(--run "${RUN_NAME}" --cam-mode "${cam_mode}")
        if [[ -n "${MOTIF_DB}" ]]; then
            MODISCO_ARGS+=(--motif-db "${MOTIF_DB}")
        fi
        PYTHON_BIN="${TRAIN_PYTHON}" "${SCRIPT_DIR}/run_modisco_report.sh" "${MODISCO_ARGS[@]}"
    done
fi

if [[ "${DO_PLOTS}" -eq 1 ]]; then
    IFS=',' read -r -a CAM_MODE_LIST <<< "${CAM_MODES}"
    for cam_mode in "${CAM_MODE_LIST[@]}"; do
        cam_mode="$(echo "${cam_mode}" | xargs)"
        [[ -z "${cam_mode}" ]] && continue
        PLOT_ARGS=(--run "${RUN_NAME}" --cam-mode "${cam_mode}")
        if [[ -n "${FOLD}" ]]; then
            PLOT_ARGS+=(--fold "${FOLD}")
        fi
        if [[ -n "${CHECKPOINT}" ]]; then
            PLOT_ARGS+=(--checkpoint "${CHECKPOINT}")
        fi
        "${TRAIN_PYTHON}" "${PLOT_EXPORT_PY}" "${PLOT_ARGS[@]}"
    done
fi
