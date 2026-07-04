#!/usr/bin/env bash

set -euo pipefail

# Logging

log() {
    local msg="$1"
    printf "[%s] %s\n" "$(date '+%F %T')" "$msg" >&2
}

warn() {
    local msg="$1"
    printf "[%s] [WARNING] %s\n" "$(date '+%F %T')" "$msg" >&2
}

fail() {
    local msg="$1"
    printf "[%s] [ERROR] %s\n" "$(date '+%F %T')" "$msg" >&2
    exit 1
}

start_step_log() {
    local log_file="$1"

    mkdir -p "$(dirname "$log_file")"

    # 同时输出到终端和日志文件
    exec > >(tee -a "$log_file") 2>&1

    log "Log file: $log_file"
}

# Dependency checks

require_cmd() {
    local cmd="$1"

    command -v "$cmd" >/dev/null 2>&1 ||
        fail "Command not found: $cmd"
}

check_dependencies() {
    local cmd

    for cmd in "$@"; do
        require_cmd "$cmd"
    done
}

# Path utilities

resolve_project_path() {
    local path="$1"

    if [[ "$path" = /* ]]; then
        printf "%s\n" "$path"
    else
        printf "%s/%s\n" "$PROJECT_DIR" "$path"
    fi
}

# File checks

check_file() {
    local file="$1"

    [[ -s "$file" ]] ||
        fail "File missing or empty: $file"
}

check_dir() {
    local dir="$1"

    [[ -d "$dir" ]] ||
        fail "Directory not found: $dir"
}

make_dir() {
    mkdir -p "$1"
}

# FASTQ checks

check_fastq() {
    local fq="$1"

    check_file "$fq"

    case "$fq" in
        *.fastq | *.fq | *.fastq.gz | *.fq.gz)
            ;;
        *)
            fail "Not a FASTQ file: $fq"
            ;;
    esac
}

# Sample table

check_sample_table() {
    local table="$1"
    local header
    local expected_header

    check_file "$table"

    header="$(head -n 1 "$table" | tr -d '\r')"
    expected_header=$'run_id\tcell_line\ttime\ttreatment\tfq'

    [[ "$header" == "$expected_header" ]] ||
        fail "Invalid sample table header.
Expected:
run_id<TAB>cell_line<TAB>time<TAB>treatment<TAB>fq
Found:
$header"
}

read_sample_table() {
    local table="$1"

    check_sample_table "$table"
    tail -n +2 "$table"
}

# Command execution

run_cmd() {
    local log_file="$1"
    shift

    mkdir -p "$(dirname "$log_file")"

    log "Running: $*"

    if ! "$@" >"$log_file" 2>&1; then
        tail -n 30 "$log_file" >&2 || true
        fail "Command failed. See log: $log_file"
    fi
}

# HISAT2 index check

check_hisat2_index() {
    local prefix="$1"
    local extension=""
    local i

    if [[ -s "${prefix}.1.ht2" ]]; then
        extension="ht2"
    elif [[ -s "${prefix}.1.ht2l" ]]; then
        extension="ht2l"
    else
        fail "HISAT2 index not found: ${prefix}.[1-8].ht2"
    fi

    for i in {1..8}; do
        check_file "${prefix}.${i}.${extension}"
    done
}

# BAM and count checks

check_bam() {
    local bam="$1"

    check_file "$bam"

    samtools quickcheck "$bam" ||
        fail "Invalid BAM file: $bam"
}

check_counts() {
    local counts="$1"
    local n

    check_file "$counts"

    n="$(wc -l <"$counts")"

    [[ "$n" -gt 1 ]] ||
        fail "Empty count table: $counts"
}