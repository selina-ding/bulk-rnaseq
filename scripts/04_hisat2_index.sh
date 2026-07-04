#!/usr/bin/env bash

set -euo pipefail

source config/config.sh
source scripts/lib.sh

mkdir -p "${REF_DIR}/hisat2_index"

hisat2-build \
    -p "$THREADS" \
    data/reference/GRCh38.primary_assembly.genome.fa \
    data/reference/hisat2_index/GRCh38 \
    > logs/hisat2-build.log 2>&1

log "HISAT2 index finished"