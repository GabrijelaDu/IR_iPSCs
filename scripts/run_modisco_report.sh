#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PIPELINE_PY="${SCRIPT_DIR}/hl_revision_pipeline.py"
PARNET_PYTHON="/lustre/groups/crna01/projects/collabs/shared/cenvs/parnet_clean/bin/python"
PYTHON_BIN="${PYTHON_BIN:-${PARNET_PYTHON}}"

usage() {
    cat <<'EOF'
Usage:
  run_modisco_report.sh --run RUN --motif-db PATH [--side left,right] [--window 256] [--size 31] [--max-seqlets 1000] [--cam-mode MODE]

Notes:
  - Use `parnet_clean` for validation and export steps; this script only uses `tfmodisco` for motif/report execution.
  - The tf-modisco inputs must already exist under the run's interpretation directory.
  - The motif database path is required because the older hardcoded repo path was stale.
EOF
}

RUN_NAME=""
MOTIF_DB=""
WINDOW_SIZE=256
MOTIF_SIZE=31
MAX_SEQLETS=1000
SIDES="left,right"
RUN_REPORT=1
CAM_MODE="final_logit_linearized"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run)
            RUN_NAME="$2"
            shift 2
            ;;
        --motif-db)
            MOTIF_DB="$2"
            shift 2
            ;;
        --side|--sides)
            SIDES="$2"
            shift 2
            ;;
        --window)
            WINDOW_SIZE="$2"
            shift 2
            ;;
        --size)
            MOTIF_SIZE="$2"
            shift 2
            ;;
        --max-seqlets|--n-seqs)
            MAX_SEQLETS="$2"
            shift 2
            ;;
        --skip-report)
            RUN_REPORT=0
            shift
            ;;
        --cam-mode)
            CAM_MODE="$2"
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

if [[ "${RUN_REPORT}" -eq 1 && -z "${MOTIF_DB}" ]]; then
    echo "A motif database is required for the report step." >&2
    exit 1
fi

eval "$("${PYTHON_BIN}" "${PIPELINE_PY}" resolve-run --run "${RUN_NAME}" --format shell)"

if [[ ! -x "${MODISCO_PYTHON}" ]]; then
    echo "tf-modisco python not found: ${MODISCO_PYTHON}" >&2
    exit 1
fi
if [[ ! -f "${MODISCO_BIN}" ]]; then
    echo "tf-modisco executable not found: ${MODISCO_BIN}" >&2
    exit 1
fi
if [[ "${RUN_REPORT}" -eq 1 && ! -f "${MOTIF_DB}" ]]; then
    echo "Motif database not found: ${MOTIF_DB}" >&2
    exit 1
fi

MODISCO_BINDIR="$(dirname "${MODISCO_PYTHON}")"

case "${CAM_MODE}" in
    final_logit_linearized)
        MODE_PREFIX="${MODISCO_PREFIX}"
        ;;
    branch_signed)
        MODE_PREFIX="${MODISCO_PREFIX}_${CAM_MODE}"
        ;;
    *)
        echo "Unsupported CAM mode: ${CAM_MODE}" >&2
        exit 1
        ;;
esac

prepare_inputs_for_modisco() {
    local cam_file="$1"
    local seq_file="$2"
    local output_dir="$3"
    local side_name="$4"

    "${PYTHON_BIN}" - "${cam_file}" "${seq_file}" "${output_dir}" "${side_name}" <<'PY'
import sys
from pathlib import Path

import numpy as np

cam_path = Path(sys.argv[1])
seq_path = Path(sys.argv[2])
out_dir = Path(sys.argv[3])
side_name = sys.argv[4]
arr = np.load(cam_path)["arr_0"].astype(np.float32)
seq = np.load(seq_path)["arr_0"].astype(np.float32)
reasons = []

if arr.shape != seq.shape:
    raise ValueError(f"Sequence/CAM shape mismatch for {side_name}: {seq.shape} vs {arr.shape}")

if not np.isfinite(arr).all():
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    reasons.append("replaced non-finite CAM values with 0")
if not np.isfinite(seq).all():
    seq = np.nan_to_num(seq, nan=0.0, posinf=0.0, neginf=0.0)
    reasons.append("replaced non-finite sequence values with 0")

seq_nlc = np.moveaxis(seq, 1, 2)
arr_nlc = np.moveaxis(arr, 1, 2)
inactive_mask = seq_nlc.sum(axis=2) <= 0
if inactive_mask.any():
    seq_nlc[inactive_mask] = 0.25
    arr_nlc[inactive_mask] = 0.0
    reasons.append(f"replaced {int(inactive_mask.sum())} inactive sequence rows with uniform background")
seq = np.moveaxis(seq_nlc, 2, 1)
arr = np.moveaxis(arr_nlc, 2, 1)

track = arr.sum(axis=1)
active_mask = seq.sum(axis=1) > 0
active_values = track[active_mask]

if active_values.size == 0:
    print(f"Skipping {side_name}: attribution file contains no active positions.", file=sys.stderr)
    print("")
    raise SystemExit(0)

has_pos = bool((active_values > 0).any())
has_neg = bool((active_values < 0).any())

adjusted_track = track.copy()
if not (has_pos and has_neg):
    if not has_pos and not has_neg:
        print(f"Skipping {side_name}: attribution file contains only zeros.", file=sys.stderr)
        print("")
        raise SystemExit(0)

    abs_active = np.abs(active_values)
    base_scale = float(np.quantile(abs_active, 0.05)) if active_values.size > 20 else float(np.median(abs_active))
    base_scale = max(base_scale, 1e-6)

    if has_pos and not has_neg:
        adjusted_track[active_mask] = adjusted_track[active_mask] - np.float32(base_scale)
        if not (adjusted_track[active_mask] < 0).any():
            low_q = float(np.quantile(active_values, 0.25))
            adjusted_track[active_mask] = adjusted_track[active_mask] - np.float32(low_q + base_scale)
        reasons.append(f"subtracted a low positive baseline ({base_scale:.3g}) from active positions")
    else:
        adjusted_track[active_mask] = adjusted_track[active_mask] + np.float32(base_scale)
        if not (adjusted_track[active_mask] > 0).any():
            high_q = float(np.quantile(active_values, 0.75))
            adjusted_track[active_mask] = adjusted_track[active_mask] - np.float32(high_q - base_scale)
        reasons.append(f"added a low negative baseline ({base_scale:.3g}) to active positions")

    if not ((adjusted_track[active_mask] > 0).any() and (adjusted_track[active_mask] < 0).any()):
        jitter = np.float32(max(base_scale * 0.05, 1e-6))
        cutoff = np.quantile(np.abs(adjusted_track[active_mask]), 0.25)
        low_signal_mask = active_mask & (np.abs(adjusted_track) <= cutoff)
        low_rows, low_cols = np.where(low_signal_mask)
        if low_rows.size:
            signs = np.where(np.arange(low_rows.size) % 2 == 0, -1.0, 1.0).astype(np.float32)
            adjusted_track[low_rows, low_cols] = signs * jitter
            reasons.append(f"added deterministic low-signal jitter ({jitter:.3g})")

    if not ((adjusted_track[active_mask] > 0).any() and (adjusted_track[active_mask] < 0).any()):
        raise ValueError(f"Could not prepare one-sided CAM track for tf-modisco on {side_name}.")

adjusted = np.zeros_like(arr, dtype=np.float32)
carrier = seq.argmax(axis=1)
rows, cols = np.where(active_mask)
adjusted[rows, carrier[rows, cols], cols] = adjusted_track[rows, cols]

cam_out_path = out_dir / f"{cam_path.stem}.modisco_ready.npz"
seq_out_path = out_dir / f"{seq_path.stem}.modisco_ready.npz"
np.savez_compressed(cam_out_path, arr_0=adjusted)
np.savez_compressed(seq_out_path, arr_0=seq)
if reasons:
    print(f"Prepared {side_name} inputs for tf-modisco: {'; '.join(reasons)}.", file=sys.stderr)
print(f"{cam_out_path}\t{seq_out_path}")
PY
}

IFS=',' read -r -a SIDE_LIST <<< "${SIDES}"

echo "Running tf-modisco for ${RUN_NAME}"
echo "Resolver python: ${PYTHON_BIN}"
echo "TF-MoDISco python: ${MODISCO_PYTHON}"
echo "TF-MoDISco bin dir: ${MODISCO_BINDIR}"
echo "Interpretation directory: ${INTERPRETATION_DIR}"
echo "Window=${WINDOW_SIZE} Size=${MOTIF_SIZE} MaxSeqlets=${MAX_SEQLETS}"
echo "CAM mode: ${CAM_MODE}"

for side in "${SIDE_LIST[@]}"; do
    side="$(echo "${side}" | xargs)"
    seq_file="${INTERPRETATION_DIR}/${MODE_PREFIX}_seq_${side}.npz"
    cam_file="${INTERPRETATION_DIR}/${MODE_PREFIX}_cam_${side}.npz"
    output_dir="${INTERPRETATION_DIR}/${MODE_PREFIX}_${side}_modisco"
    output_h5="${output_dir}/${MODE_PREFIX}_${side}_modisco.h5"

    if [[ ! -f "${seq_file}" ]]; then
        echo "Missing sequence file: ${seq_file}" >&2
        exit 1
    fi
    if [[ ! -f "${cam_file}" ]]; then
        echo "Missing CAM file: ${cam_file}" >&2
        exit 1
    fi

    mkdir -p "${output_dir}"

    IFS=$'\t' read -r prepared_cam_file prepared_seq_file <<< "$(prepare_inputs_for_modisco "${cam_file}" "${seq_file}" "${output_dir}" "${side}")"
    if [[ -z "${prepared_cam_file}" || -z "${prepared_seq_file}" ]]; then
        continue
    fi

    PATH="${MODISCO_BINDIR}:${PATH}" "${MODISCO_PYTHON}" "${MODISCO_BIN}" motifs \
        -s "${prepared_seq_file}" \
        -a "${prepared_cam_file}" \
        -n "${MAX_SEQLETS}" \
        --window "${WINDOW_SIZE}" \
        --size "${MOTIF_SIZE}" \
        -o "${output_h5}"

    if [[ ! -f "${output_h5}" ]]; then
        echo "tf-modisco motifs step did not create ${output_h5}" >&2
        exit 1
    fi

    if [[ "${RUN_REPORT}" -eq 1 ]]; then
        report_dir="${output_dir}/report"
        mkdir -p "${report_dir}"
        PATH="${MODISCO_BINDIR}:${PATH}" "${MODISCO_PYTHON}" "${MODISCO_BIN}" report \
            -i "${output_h5}" \
            -o "${report_dir}" \
            -s "${report_dir}" \
            -m "${MOTIF_DB}"
        echo "Report written to ${report_dir}"
    else
        echo "Motif discovery written to ${output_h5}"
    fi
done
