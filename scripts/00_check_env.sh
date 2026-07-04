#!/usr/bin/env bash
source "$(dirname "$0")/../config/config.sh"
source "$(dirname "$0")/lib.sh"

for cmd in \
    fasterq-dump \
    fastqc \
    multiqc \
    fastp \
    hisat2 \
    samtools \
    featureCounts \
    Rscript \
    snakemake
do
    require_cmd "$cmd"
done

log "All required commands are available."
log "Project directory: $PROJECT_DIR"