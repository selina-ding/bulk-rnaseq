#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/lib.sh"

[[ $# -ge 1 ]] ||
    fail "Usage: bash scripts/06_featurecounts.sh <sample_table.tsv> [run_name]"

SAMPLE_TABLE="$(resolve_project_path "$1")"
RUN_NAME="${2:-$(basename "$SAMPLE_TABLE" .tsv)}"

BAM_INPUT="${BAM_DIR}/${RUN_NAME}"
COUNTS_OUT="${COUNTS_DIR}/${RUN_NAME}"
STEP_LOG_DIR="${LOG_DIR}/06_featurecounts/${RUN_NAME}"

mkdir -p "$COUNTS_OUT" "$STEP_LOG_DIR"

start_step_log "${STEP_LOG_DIR}/06_featurecounts.log"

check_dependencies featureCounts samtools
check_sample_table "$SAMPLE_TABLE"
check_file "$GTF_FILE"

BAM_FILES=()

while IFS=$'\t' read -r run_id cell_line time treatment fq
do
    [[ -z "$run_id" ]] && continue

    bam_file="${BAM_INPUT}/${run_id}.sorted.bam"

    check_bam "$bam_file"

    BAM_FILES+=("$bam_file")

done < <(read_sample_table "$SAMPLE_TABLE")

[[ "${#BAM_FILES[@]}" -gt 0 ]] ||
    fail "No BAM files were found for run: $RUN_NAME"

COUNT_FILE="${COUNTS_OUT}/gene_counts.tsv"
COMMAND_LOG="${STEP_LOG_DIR}/featurecounts-command.log"

log "Running featureCounts with ${#BAM_FILES[@]} BAM files"

if ! featureCounts \
    -T "$THREADS" \
    -a "$GTF_FILE" \
    -t exon \
    -g gene_id \
    -s "$STRANDEDNESS" \
    -o "$COUNT_FILE" \
    "${BAM_FILES[@]}" \
    >"$COMMAND_LOG" 2>&1
then
    tail -n 30 "$COMMAND_LOG" >&2 || true
    fail "featureCounts failed: $RUN_NAME"
fi

check_counts "$COUNT_FILE"
check_file "${COUNT_FILE}.summary"

log "featureCounts finished: $RUN_NAME"
log "Count table: $COUNT_FILE"