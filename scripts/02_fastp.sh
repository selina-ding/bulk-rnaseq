#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/lib.sh"

[[ $# -ge 1 ]] ||
    fail "Usage: bash scripts/02_fastp.sh <sample_table.tsv> [run_name]"

SAMPLE_TABLE="$(resolve_project_path "$1")"
RUN_NAME="${2:-$(basename "$SAMPLE_TABLE" .tsv)}"

TRIM_OUT="${TRIM_DIR}/${RUN_NAME}"
FASTP_QC_OUT="${QC_DIR}/${RUN_NAME}/fastp"
STEP_LOG_DIR="${LOG_DIR}/02_fastp/${RUN_NAME}"

mkdir -p "$TRIM_OUT" "$FASTP_QC_OUT" "$STEP_LOG_DIR"

start_step_log "${STEP_LOG_DIR}/02_fastp.log"

check_dependencies fastp
check_sample_table "$SAMPLE_TABLE"

log "Starting fastp: $RUN_NAME"

while IFS=$'\t' read -r run_id cell_line time treatment fq
do
    [[ -z "$run_id" ]] && continue

    fq_path="$(resolve_project_path "$fq")"

    check_fastq "$fq_path"

    out_fq="${TRIM_OUT}/${run_id}.trimmed.fastq.gz"
    json="${FASTP_QC_OUT}/${run_id}.fastp.json"
    html="${FASTP_QC_OUT}/${run_id}.fastp.html"
    sample_log="${STEP_LOG_DIR}/${run_id}.log"

    log "Processing: $run_id"

    if ! fastp \
        --in1 "$fq_path" \
        --out1 "$out_fq" \
        --thread "$THREADS" \
        --json "$json" \
        --html "$html" \
        >"$sample_log" 2>&1
    then
        tail -n 30 "$sample_log" >&2 || true
        fail "fastp failed: $run_id"
    fi

    check_fastq "$out_fq"

done < <(read_sample_table "$SAMPLE_TABLE")

log "fastp finished: $RUN_NAME"