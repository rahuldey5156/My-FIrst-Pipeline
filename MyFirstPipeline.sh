#!/usr/bin/env bash

# =============================================================================
# MyFirstPipeline.sh
# RNAseq analysis pipeline for Trypanosoma congolense RNAi experiment
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR="/localdisk/data/BPSM/MyFirstPipeline"
FASTQ_DIR="${BASE_DIR}/fastq"
SAMPLE_SHEET="${BASE_DIR}/Tco2.fqfiles"
GENOME_DIR="${BASE_DIR}/Tcongo_genome"
BED_FILE="${BASE_DIR}/TriTrypDB-46_TcongolenseIL3000_2019.bed"

OUT_DIR="${BASE_DIR}/pipeline_output"
QC_DIR="${OUT_DIR}/fastqc"
INDEX_DIR="${OUT_DIR}/bowtie2_index"
BAM_DIR="${OUT_DIR}/bam"
COUNTS_DIR="${OUT_DIR}/counts"
MEANS_DIR="${OUT_DIR}/means"
FC_DIR="${OUT_DIR}/fold_changes"
LOG_DIR="${OUT_DIR}/logs"

INDEX_PREFIX="${INDEX_DIR}/TcongolenseIL3000"
GENOME_CAT="${INDEX_DIR}/TcongolenseIL3000_genome.fa"

THREADS=4

# =============================================================================
# LOGGING
# =============================================================================
# All messages go to stdout and a timestamped log file in OUT_DIR.

LOG_FILE="${OUT_DIR}/pipeline_run_$(date +%Y%m%d_%H%M%S).log"

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg"
    if [[ "${LOG_READY:-0}" == "1" ]]; then
        echo "$msg" >> "$LOG_FILE"
    fi
}

# =============================================================================
# CREATE OUTPUT DIRECTORIES
# =============================================================================

mkdir -p "$QC_DIR" "$INDEX_DIR" "$BAM_DIR" "$COUNTS_DIR" "$MEANS_DIR" "$FC_DIR" "$LOG_DIR"
LOG_READY=1
log "Output directories created under ${OUT_DIR}"
