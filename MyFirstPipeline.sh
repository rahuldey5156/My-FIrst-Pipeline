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

# =============================================================================
# STEP 1 – SANITY CHECKS
# =============================================================================
# Verify all required tools exist in PATH, input files/dirs are present,
# the sample sheet has the right format, all fastq files exist, and that
# at least one genome FASTA is available.

log "========== STEP 1: Sanity checks =========="

# -- 1a. Required tools --
REQUIRED_TOOLS=(fastqc bowtie2 bowtie2-build samtools bedtools awk sort)
for tool in "${REQUIRED_TOOLS[@]}"; do
    if ! command -v "$tool" &>/dev/null; then
        log "ERROR: Required tool '$tool' not found in PATH. Aborting."
        exit 1
    fi
    log "  Tool found: $tool"
done

# -- 1b. Required paths --
for path in "$FASTQ_DIR" "$GENOME_DIR" "$SAMPLE_SHEET" "$BED_FILE"; do
    if [[ ! -e "$path" ]]; then
        log "ERROR: Required path does not exist: $path"
        exit 1
    fi
done
log "  All required input paths verified."

# -- 1c. Sample sheet non-empty and has at least 4 tab-separated columns --
if [[ ! -s "$SAMPLE_SHEET" ]]; then
    log "ERROR: Sample sheet is empty: $SAMPLE_SHEET"
    exit 1
fi

awk 'NF < 4 { print NR": "$0; err=1 } END { exit err+0 }' "$SAMPLE_SHEET" || {
    log "ERROR: One or more lines in the sample sheet have fewer than 4 fields."
    exit 1
}
log "  Sample sheet format OK."

# -- 1d. All fastq files referenced in sample sheet actually exist --
while IFS=$'\t' read -r sample r1 r2 group _rest; do
    [[ "$sample" =~ ^#.*$ ]] && continue
    for fq in "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2"; do
        if [[ ! -f "$fq" ]]; then
            log "ERROR: fastq file not found: $fq  (sample: '$sample')"
            exit 1
        fi
    done
done < "$SAMPLE_SHEET"
log "  All fastq files present."

# -- 1e. At least one genome FASTA file exists --
genome_files=( $(ls "$GENOME_DIR"/*.fa "$GENOME_DIR"/*.fasta 2>/dev/null || true) )
if [[ ${#genome_files[@]} -eq 0 ]]; then
    log "ERROR: No .fa or .fasta files found in: $GENOME_DIR"
    exit 1
fi
log "  Genome files found: ${#genome_files[@]} file(s)."
log "  All sanity checks passed."

# =============================================================================
# STEP 2 – QC #1: FastQC on all raw fastq.gz files
# =============================================================================
# Runs FastQC on every read file listed in the sample sheet.
# Skips files whose output zip already exists (safe to re-run).

log "========== STEP 2: FastQC =========="

while IFS=$'\t' read -r sample r1 r2 group _rest; do
    [[ "$sample" =~ ^#.*$ ]] && continue
    for fq in "$r1" "$r2"; do
        fq_path="${FASTQ_DIR}/${fq}"
        base="${fq%.fq.gz}"; base="${base%.fastq.gz}"
        zip_out="${QC_DIR}/${base}_fastqc.zip"

        if [[ -f "$zip_out" ]]; then
            log "  FastQC output exists, skipping: $fq"
        else
            log "  Running FastQC on: $fq"
            fastqc \
                --outdir "$QC_DIR" \
                --threads "$THREADS" \
                --quiet \
                "$fq_path" \
                >> "${LOG_DIR}/fastqc.log" 2>&1
            log "  FastQC done: $fq"
        fi
    done
done < "$SAMPLE_SHEET"

log "  FastQC reports written to: $QC_DIR"

# =============================================================================
# STEP 3 – QC #2: Parse FastQC summaries; report PASS/WARN/FAIL per sample
# =============================================================================
# Each FastQC zip contains summary.txt. We extract it, collate all results
# into a single QC_summary.txt, and warn loudly if any module FAILs.

log "========== STEP 3: Summarise FastQC results =========="

QC_SUMMARY="${OUT_DIR}/QC_summary.txt"
printf "%-50s %-40s %s\n" "Sample" "Module" "Result" > "$QC_SUMMARY"
printf "%s\n" "$(printf '=%.0s' {1..100})"            >> "$QC_SUMMARY"

FAIL_COUNT=0

for zip in "$QC_DIR"/*_fastqc.zip; do
    [[ -f "$zip" ]] || continue
    tmpdir=$(mktemp -d)
    unzip -q "$zip" "*/summary.txt" -d "$tmpdir" 2>/dev/null || {
        log "  WARNING: Could not extract summary.txt from $zip"
        rm -rf "$tmpdir"; continue
    }
    summary_file=$(find "$tmpdir" -name "summary.txt" | head -1)

    while IFS=$'\t' read -r result module filename; do
        printf "%-50s %-40s %s\n" "$filename" "$module" "$result" >> "$QC_SUMMARY"
        if [[ "$result" == "FAIL" ]]; then
            log "  QC FAIL – $filename : $module"
            (( FAIL_COUNT++ )) || true
        fi
    done < "$summary_file"
    rm -rf "$tmpdir"
done

log "  QC summary written to: $QC_SUMMARY"
if [[ "$FAIL_COUNT" -gt 0 ]]; then
    log "  WARNING: $FAIL_COUNT module(s) FAILED FastQC. Review $QC_SUMMARY."
else
    log "  No hard FAILs in FastQC output."
fi

# =============================================================================
# STEP 4 – Build bowtie2 genome index (once only)
# =============================================================================
# The genome is split across multiple per-chromosome FASTA files; we first
# concatenate them into one file, then run bowtie2-build.
# A sentinel file (.index_built) prevents needless re-indexing on re-runs.

log "========== STEP 4: Build bowtie2 genome index =========="

INDEX_STAMP="${INDEX_DIR}/.index_built"

if [[ -f "$INDEX_STAMP" ]]; then
    log "  Index already built (delete ${INDEX_STAMP} to force rebuild). Skipping."
else
    log "  Concatenating genome FASTAs into: $GENOME_CAT"
    cat "$GENOME_DIR"/*.fa "$GENOME_DIR"/*.fasta > "$GENOME_CAT" 2>/dev/null || \
    cat "$GENOME_DIR"/*.fa  > "$GENOME_CAT" 2>/dev/null || \
    cat "$GENOME_DIR"/*.fasta > "$GENOME_CAT"

    [[ -s "$GENOME_CAT" ]] || { log "ERROR: Concatenated genome FASTA is empty."; exit 1; }

    log "  Running bowtie2-build …"
    bowtie2-build \
        --threads "$THREADS" \
        "$GENOME_CAT" \
        "$INDEX_PREFIX" \
        >> "${LOG_DIR}/bowtie2_build.log" 2>&1

    touch "$INDEX_STAMP"
    log "  Index built: ${INDEX_PREFIX}.*"
fi

# =============================================================================
# STEP 5 – Align reads with bowtie2; convert to sorted indexed BAM
# =============================================================================
# For each sample:
#   bowtie2         → paired-end alignment (SAM piped to samtools)
#   samtools view   → SAM → BAM, discard unmapped reads (-F 4)
#   samtools sort   → coordinate-sort the BAM
#   samtools index  → create .bai index required by bedtools
#
# bowtie2 flags:
#   --no-unal        suppress unaligned reads from output (smaller files)
#   --no-mixed       require both mates to align
#   --no-discordant  require concordant pairs only (correct orientation/distance)

log "========== STEP 5: Align reads with bowtie2 =========="

while IFS=$'\t' read -r sample r1 r2 group _rest; do
    [[ "$sample" =~ ^#.*$ ]] && continue

    bam_out="${BAM_DIR}/${sample}.sorted.bam"

    if [[ -f "${bam_out}.bai" ]]; then
        log "  BAM already exists for ${sample}, skipping."
        continue
    fi

    log "  Aligning: ${sample}  (group: ${group})"

    bowtie2 \
        -x "$INDEX_PREFIX" \
        -1 "${FASTQ_DIR}/${r1}" \
        -2 "${FASTQ_DIR}/${r2}" \
        --no-unal \
        --no-mixed \
        --no-discordant \
        -p "$THREADS" \
        2>> "${LOG_DIR}/${sample}_bowtie2.log" \
    | samtools view -bS -F 4 - \
    | samtools sort -@ "$THREADS" -o "$bam_out" -

    samtools index "$bam_out"
    log "  Done: ${bam_out}"
done < "$SAMPLE_SHEET"

log "  All alignments complete. BAMs in: $BAM_DIR"
