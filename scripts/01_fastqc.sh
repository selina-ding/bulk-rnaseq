#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/lib.sh"

[[ $# -ge 1 ]] ||
    fail "Usage: bash scripts/01_fastqc.sh <sample_table.tsv> [run_name]"

SAMPLE_TABLE="$(resolve_project_path "$1")"
RUN_NAME="${2:-$(basename "$SAMPLE_TABLE" .tsv)}"

FASTQC_OUT="${QC_DIR}/${RUN_NAME}/fastqc"
STEP_LOG_DIR="${LOG_DIR}/01_fastqc/${RUN_NAME}"

mkdir -p "$FASTQC_OUT" "$STEP_LOG_DIR"

start_step_log "${STEP_LOG_DIR}/01_fastqc.log"

check_dependencies fastqc
check_sample_table "$SAMPLE_TABLE"

log "Starting FastQC: $RUN_NAME"

while IFS=$'\t' read -r run_id cell_line time treatment fq
do
    [[ -z "$run_id" ]] && continue

    fq_path="$(resolve_project_path "$fq")"

    check_fastq "$fq_path"

    log "Processing: $run_id"

    if ! fastqc \
        "$fq_path" \
        --threads "$THREADS" \
        --outdir "$FASTQC_OUT" \
        >"${STEP_LOG_DIR}/${run_id}.log" 2>&1
    then
        tail -n 30 "${STEP_LOG_DIR}/${run_id}.log" >&2 || true
        fail "FastQC failed: $run_id"
    fi

done < <(read_sample_table "$SAMPLE_TABLE")

log "FastQC finished: $RUN_NAME"