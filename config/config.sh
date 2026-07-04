#!/usr/bin/env bash

########################################
# Project root
########################################

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

########################################
# Computing resources
########################################

THREADS="${THREADS:-8}"
MEM_GB="${MEM_GB:-32}"

# featureCounts:
# 0 = unstranded
# 1 = stranded
# 2 = reversely stranded
STRANDEDNESS="${STRANDEDNESS:-0}"

########################################
# Data directories
########################################

DATA_DIR="${PROJECT_DIR}/data"
RAW_DIR="${DATA_DIR}/raw_fastq"
TRIM_DIR="${DATA_DIR}/trimmed_fastq"
REF_DIR="${DATA_DIR}/reference"

########################################
# Reference files
########################################

GENOME_FA="${REF_DIR}/GRCh38.primary_assembly.genome.fa"
GTF_FILE="${REF_DIR}/gencode.v44.annotation.gtf"

HISAT2_INDEX_DIR="${REF_DIR}/hisat2_index"
HISAT2_INDEX="${HISAT2_INDEX_DIR}/GRCh38"

########################################
# Result directories
########################################

RESULTS_DIR="${PROJECT_DIR}/results"

QC_DIR="${RESULTS_DIR}/qc"
BAM_DIR="${RESULTS_DIR}/bam"
COUNTS_DIR="${RESULTS_DIR}/counts"
DESEQ2_DIR="${RESULTS_DIR}/deseq2"
ENRICH_DIR="${RESULTS_DIR}/enrichment"
FIGURES_DIR="${RESULTS_DIR}/figures"
GSEA_DIR="${RESULTS_DIR}/gsea"
REPORT_DIR="${RESULTS_DIR}/reports"

########################################
# Log directory
########################################

LOG_DIR="${PROJECT_DIR}/logs"

########################################
# Create base directories
########################################

mkdir -p \
    "$RAW_DIR" \
    "$TRIM_DIR" \
    "$REF_DIR" \
    "$HISAT2_INDEX_DIR" \
    "$QC_DIR" \
    "$BAM_DIR" \
    "$COUNTS_DIR" \
    "$DESEQ2_DIR" \
    "$ENRICH_DIR" \
    "$FIGURES_DIR" \
    "$GSEA_DIR" \
    "$REPORT_DIR" \
    "$LOG_DIR"