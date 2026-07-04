#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/lib.sh"

[[ $# -ge 1 ]] ||
    fail "Usage: bash scripts/03_multiqc.sh <sample_table.tsv> [run_name]"

SAMPLE_TABLE="$(resolve_project_path "$1")"
RUN_NAME="${2:-$(basename "$SAMPLE_TABLE" .tsv)}"

QC_INPUT="${QC_DIR}/${RUN_NAME}"
MULTIQC_OUT="${QC_INPUT}/multiqc"
STEP_LOG_DIR="${LOG_DIR}/03_multiqc/${RUN_NAME}"

mkdir -p "$MULTIQC_OUT" "$STEP_LOG_DIR"

start_step_log "${STEP_LOG_DIR}/03_multiqc.log"

check_dependencies multiqc
check_sample_table "$SAMPLE_TABLE"
check_dir "$QC_INPUT"

log "Starting MultiQC: $RUN_NAME"

if ! multiqc \
    "$QC_INPUT" \
    --outdir "$MULTIQC_OUT" \
    --filename "multiqc_report.html" \
    --force \
    >"${STEP_LOG_DIR}/multiqc-command.log" 2>&1
then
    tail -n 30 "${STEP_LOG_DIR}/multiqc-command.log" >&2 || true
    fail "MultiQC failed: $RUN_NAME"
fi

check_file "${MULTIQC_OUT}/multiqc_report.html"

log "MultiQC finished: $RUN_NAME"