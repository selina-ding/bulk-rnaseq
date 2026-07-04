#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "${SCRIPT_DIR}/../config/config.sh"
source "${SCRIPT_DIR}/lib.sh"

[[ $# -ge 1 ]] ||
    fail "Usage: bash scripts/05_hisat2_align.sh <sample_table.tsv> [run_name]"

SAMPLE_TABLE="$(resolve_project_path "$1")"
RUN_NAME="${2:-$(basename "$SAMPLE_TABLE" .tsv)}"

TRIM_INPUT="${TRIM_DIR}/${RUN_NAME}"
BAM_OUT="${BAM_DIR}/${RUN_NAME}"
STEP_LOG_DIR="${LOG_DIR}/05_hisat2_align/${RUN_NAME}"

# Conservative defaults for a computer with limited memory.
# You can override these values when running the script.
HISAT2_THREADS="${HISAT2_THREADS:-2}"
SORT_THREADS="${SORT_THREADS:-1}"
SORT_MEMORY="${SORT_MEMORY:-768M}"

mkdir -p "$BAM_OUT" "$STEP_LOG_DIR"

start_step_log "${STEP_LOG_DIR}/05_hisat2_align.log"

check_dependencies hisat2 samtools
check_sample_table "$SAMPLE_TABLE"
check_hisat2_index "$HISAT2_INDEX"

log "Starting HISAT2 alignment: $RUN_NAME"
log "HISAT2 threads: $HISAT2_THREADS"
log "samtools sort threads: $SORT_THREADS"
log "samtools sort memory per thread: $SORT_MEMORY"

while IFS=$'\t' read -r run_id cell_line time treatment fq
do
    [[ -z "$run_id" ]] && continue

    trimmed_fq="${TRIM_INPUT}/${run_id}.trimmed.fastq.gz"

    bam_file="${BAM_OUT}/${run_id}.sorted.bam"
    temp_bam="${BAM_OUT}/${run_id}.sorted.tmp.bam"

    sample_log="${STEP_LOG_DIR}/${run_id}.log"
    summary_file="${STEP_LOG_DIR}/${run_id}.summary.txt"

    check_fastq "$trimmed_fq"

    log "Aligning: $run_id"

    # Remove files left by an earlier failed run.
    rm -f "$temp_bam"
    rm -f "${temp_bam}.bai"

    if ! hisat2 \
        -x "$HISAT2_INDEX" \
        -U "$trimmed_fq" \
        -p "$HISAT2_THREADS" \
        --summary-file "$summary_file" \
        2>"$sample_log" \
        | samtools sort \
            -@ "$SORT_THREADS" \
            -m "$SORT_MEMORY" \
            -T "${BAM_OUT}/${run_id}.sort_tmp" \
            -o "$temp_bam" \
            - \
            2>>"$sample_log"
    then
        rm -f "$temp_bam"
        rm -f "${BAM_OUT}/${run_id}.sort_tmp"*
        tail -n 30 "$sample_log" >&2 || true
        fail "HISAT2 alignment failed: $run_id"
    fi

    # Check the temporary BAM before replacing the final BAM.
    if ! samtools quickcheck "$temp_bam"
    then
        rm -f "$temp_bam"
        fail "Temporary BAM file is incomplete or corrupted: $run_id"
    fi

    mv "$temp_bam" "$bam_file"

    if ! samtools index "$bam_file" >>"$sample_log" 2>&1
    then
        fail "BAM indexing failed: $run_id"
    fi

    check_bam "$bam_file"
    check_file "${bam_file}.bai"

    log "Alignment completed: $run_id"

done < <(read_sample_table "$SAMPLE_TABLE")

log "HISAT2 alignment finished: $RUN_NAME"