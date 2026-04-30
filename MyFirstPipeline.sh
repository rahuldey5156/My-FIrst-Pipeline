#!/usr/bin/bash

# =============================================================================
# MyFirstPipeline.sh
# =============================================================================
# RNAseq analysis pipeline for Trypanosoma congolense RNAi knock-down experiment
# (IL3000 laboratory strain).
#
# Sample sheet columns (Tco2.fqfiles, tab-delimited):
#   $1 SampleName  $2 SampleType  $3 Replicate  $4 Time
#   $5 Treatment   $6 End1 (R1)   $7 End2 (R2)
#
# Group labels are built as: SampleType_Treatment_TTime
#   e.g. Clone1_induced_T24, WT_uninduced_T0
#
# Pipeline steps:
#   1. Sanity checks
#   2. FastQC on all raw paired-end fastq.gz files
#   3. Parse FastQC summaries → PASS/WARN/FAIL report
#   4. Build bowtie2 genome index (once only)
#   5. Align reads with bowtie2 → sorted indexed BAM
#   6. Count reads per gene with bedtools coverage
#   7. Per-group mean counts with gene descriptions
#   8. Pairwise fold-change calculations
#
# Usage:
#   bash MyFirstPipeline.sh
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

BASE_DIR="/home/s2793337/BPSM/ICA1/"
FASTQ_DIR="${BASE_DIR}/fastq"
SAMPLE_SHEET="${BASE_DIR}/fastq/Tco2.fqfiles"
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
# All messages go to stdout AND a timestamped master log file.
# LOG_READY flag is set after output dirs are created so the log
# file can be opened.

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

mkdir -p "$QC_DIR" "$INDEX_DIR" "$BAM_DIR" "$COUNTS_DIR" \
         "$MEANS_DIR" "$FC_DIR" "$LOG_DIR"
LOG_READY=1
log "Output directories created under ${OUT_DIR}"

# =============================================================================
# STEP 1 – SANITY CHECKS
# =============================================================================
# Verify tools, paths, sample sheet format, fastq files, and genome files
# are all present before running anything expensive.

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

# -- 1c. Sample sheet must be non-empty and have 7 tab-separated columns --
# The expected columns are:
#   SampleName  SampleType  Replicate  Time  Treatment  End1  End2
if [[ ! -s "$SAMPLE_SHEET" ]]; then
    log "ERROR: Sample sheet is empty: $SAMPLE_SHEET"
    exit 1
fi

awk '
    /^#/          { next }           # skip comment lines
    /^SampleName/ { next }           # skip header line
    NF < 7        { print NR": "$0; err=1 }
    END           { exit err+0 }
' "$SAMPLE_SHEET" || {
    log "ERROR: One or more data lines in the sample sheet have fewer than 7 fields."
    exit 1
}
log "  Sample sheet format OK (7 columns)."

# -- 1d. All fastq files referenced in sample sheet must exist --
# Columns 6 (End1) and 7 (End2) hold the R1 and R2 filenames.
while IFS=$'\t' read -r sample sampletype replicate time treatment r1 r2; do
    [[ "$sample" =~ ^#.*$    ]] && continue
    [[ "$sample" == "SampleName" ]] && continue
    for fq in "$FASTQ_DIR/$r1" "$FASTQ_DIR/$r2"; do
        if [[ ! -f "$fq" ]]; then
            log "ERROR: fastq file not found: $fq  (sample: '$sample')"
            exit 1
        fi
    done
done < "$SAMPLE_SHEET"
log "  All fastq files present."

# -- 1e. At least one genome FASTA (plain or gzipped) must exist --
genome_files=( $(ls "$GENOME_DIR"/*.fa    "$GENOME_DIR"/*.fasta \
                    "$GENOME_DIR"/*.fa.gz "$GENOME_DIR"/*.fasta.gz \
                    2>/dev/null || true) )
if [[ ${#genome_files[@]} -eq 0 ]]; then
    log "ERROR: No .fa / .fasta / .fa.gz / .fasta.gz files found in: $GENOME_DIR"
    exit 1
fi
log "  Genome files found: ${#genome_files[@]} file(s)."
log "  All sanity checks passed."

# =============================================================================
# STEP 2 – QC #1: FastQC on all raw fastq.gz files
# =============================================================================
# Runs FastQC on every R1 and R2 file listed in the sample sheet.
# Skips files whose output zip already exists so re-runs are safe.

log "========== STEP 2: FastQC =========="

while IFS=$'\t' read -r sample sampletype replicate time treatment r1 r2; do
    [[ "$sample" =~ ^#.*$    ]] && continue
    [[ "$sample" == "SampleName" ]] && continue

    for fq in "$r1" "$r2"; do
        fq_path="${FASTQ_DIR}/${fq}"
        # Derive expected FastQC zip name (FastQC strips .fq.gz / .fastq.gz)
        base="${fq%.fq.gz}"; base="${base%.fastq.gz}"
        zip_out="${QC_DIR}/${base}_fastqc.zip"

        if [[ -f "$zip_out" ]]; then
            log "  FastQC output exists, skipping: $fq"
        else
            log "  Running FastQC on: $fq"
            fastqc \
                --outdir  "$QC_DIR" \
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
# Each FastQC zip contains summary.txt listing module results.
# We extract and collate these into a single QC_summary.txt and warn
# loudly if any module returns FAIL.

log "========== STEP 3: Summarise FastQC results =========="

QC_SUMMARY="${OUT_DIR}/QC_summary.txt"
printf "%-50s %-40s %s\n" "Sample" "Module" "Result" >  "$QC_SUMMARY"
printf "%s\n"  "$(printf '=%.0s' {1..100})"           >> "$QC_SUMMARY"

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
# The genome may be gzipped (.fasta.gz). We decompress and concatenate all
# genome files into one FASTA, then run bowtie2-build.
# A sentinel file (.index_built) prevents needless re-indexing on re-runs.

log "========== STEP 4: Build bowtie2 genome index =========="

INDEX_STAMP="${INDEX_DIR}/.index_built"

if [[ -f "$INDEX_STAMP" ]]; then
    log "  Index already built (delete ${INDEX_STAMP} to force rebuild). Skipping."
else
    log "  Decompressing and concatenating genome FASTA(s) into: $GENOME_CAT"

    # Handle both plain and gzipped FASTA files
    for f in "$GENOME_DIR"/*.fasta.gz "$GENOME_DIR"/*.fa.gz \
              "$GENOME_DIR"/*.fasta    "$GENOME_DIR"/*.fa; do
        [[ -f "$f" ]] || continue
        case "$f" in
            *.gz) gunzip -c "$f" ;;
            *)    cat    "$f" ;;
        esac
    done > "$GENOME_CAT"

    if [[ ! -s "$GENOME_CAT" ]]; then
        log "ERROR: Concatenated genome FASTA is empty. Check $GENOME_DIR."
        exit 1
    fi

    log "  Running bowtie2-build (this may take a few minutes) …"
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
#   bowtie2        → paired-end alignment, SAM piped to samtools
#   samtools view  → SAM → BAM, discard unmapped reads (-F 4)
#   samtools sort  → coordinate-sort the BAM
#   samtools index → create .bai index required by bedtools
#
# bowtie2 flags chosen:
#   --no-unal       : suppress unaligned reads (smaller output)
#   --no-mixed      : require both mates to align
#   --no-discordant : require concordant pairs only

log "========== STEP 5: Align reads with bowtie2 =========="

while IFS=$'\t' read -r sample sampletype replicate time treatment r1 r2; do
    [[ "$sample" =~ ^#.*$    ]] && continue
    [[ "$sample" == "SampleName" ]] && continue

    # Build group label: SampleType_Treatment_TTime  e.g. Clone1_induced_T24
    group="${sampletype}_${treatment}_T${time}"
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

# =============================================================================
# STEP 6 – Count reads per gene using bedtools coverage
# =============================================================================
# bedtools coverage -counts reports how many read alignments overlap each
# gene feature in the BED file.
# Assumption: no introns, so the BED interval = full gene body.
#
# Raw output : all BED columns + count as the last column
# Clean output: gene_id (col 4) and count (last col) only — used downstream

log "========== STEP 6: Count reads per gene =========="

while IFS=$'\t' read -r sample sampletype replicate time treatment r1 r2; do
    [[ "$sample" =~ ^#.*$    ]] && continue
    [[ "$sample" == "SampleName" ]] && continue

    bam_in="${BAM_DIR}/${sample}.sorted.bam"
    counts_raw="${COUNTS_DIR}/${sample}_coverage.txt"
    counts_clean="${COUNTS_DIR}/${sample}_counts.txt"

    if [[ -f "$counts_clean" ]]; then
        log "  Counts exist for ${sample}, skipping."
        continue
    fi

    log "  Counting reads: ${sample}"

    bedtools coverage \
        -counts \
        -a "$BED_FILE" \
        -b "$bam_in" \
        >  "$counts_raw" \
        2>> "${LOG_DIR}/${sample}_bedtools.log"

    # Extract gene_id (col 4) and count (last column, NF)
    awk 'BEGIN{OFS="\t"} {print $4, $NF}' "$counts_raw" > "$counts_clean"

    log "  Counts written: $counts_clean"
done < "$SAMPLE_SHEET"

log "  All count files in: $COUNTS_DIR"

# =============================================================================
# STEP 7 – Per-group mean counts with gene descriptions
# =============================================================================
# Groups are built dynamically from the sample sheet as SampleType_Treatment_TTime,
# so no code changes are needed if new samples are added.
#
# For each group:
#   1. Collect count files for all samples in that group
#   2. AWK sums counts per gene across replicates, divides by n files
#   3. Join with gene descriptions from the BED file
#      (col 4 = gene_id, col 5 = description)
#
# Output: gene_id  mean_count  description  (one file per group)

log "========== STEP 7: Per-group mean expression levels =========="

# -- Extract gene descriptions from the BED file (done once) --
GENE_DESC="${OUT_DIR}/gene_descriptions.txt"
if [[ ! -f "$GENE_DESC" ]]; then
    awk 'BEGIN{OFS="\t"} NF>=5 {print $4, $5}' "$BED_FILE" | sort -u > "$GENE_DESC"
    log "  Gene descriptions extracted to: $GENE_DESC"
fi

# -- Detect all unique group labels from the sample sheet --
# Groups are built the same way as in the alignment/counting loops
mapfile -t GROUPS < <(
    awk '
        /^#/          { next }
        /^SampleName/ { next }
        NF >= 7       { print $2"_"$5"_T"$4 }
    ' "$SAMPLE_SHEET" | sort -u
)
log "  Groups detected: ${GROUPS[*]}"

for grp in "${GROUPS[@]}"; do

    means_out="${MEANS_DIR}/${grp}_mean_counts.txt"
    [[ -f "$means_out" ]] && { log "  Means exist for ${grp}, skipping."; continue; }

    log "  Computing means for group: ${grp}"

    # Collect count files for all samples belonging to this group
    count_files=()
    while IFS=$'\t' read -r sample sampletype replicate time treatment r1 r2; do
        [[ "$sample" =~ ^#.*$    ]] && continue
        [[ "$sample" == "SampleName" ]] && continue
        g="${sampletype}_${treatment}_T${time}"
        [[ "$g" == "$grp" ]] && count_files+=( "${COUNTS_DIR}/${sample}_counts.txt" )
    done < "$SAMPLE_SHEET"

    if [[ ${#count_files[@]} -eq 0 ]]; then
        log "  WARNING: No count files found for group ${grp}. Skipping."
        continue
    fi
    log "    ${#count_files[@]} replicate(s): ${count_files[*]}"

    # Sum counts per gene across all replicate files, then divide by n files
    awk '
    BEGIN { OFS="\t" }
    { sum[$1] += $2+0 }
    END {
        n = ARGC - 1
        for (gene in sum)
            printf "%s\t%.4f\n", gene, sum[gene]/n
    }
    ' "${count_files[@]}" | sort -k1,1 > "${MEANS_DIR}/${grp}_means_tmp.txt"

    # Join with gene descriptions; -a 1 keeps genes even if no description found
    join -t $'\t' -1 1 -2 1 -a 1 \
        "${MEANS_DIR}/${grp}_means_tmp.txt" \
        "$GENE_DESC" \
    | awk '
        BEGIN { OFS="\t"; print "gene_id","mean_count","description" }
        { print $1, $2, (NF>=3 ? $3 : "no_description") }
    ' > "$means_out"

    rm -f "${MEANS_DIR}/${grp}_means_tmp.txt"
    log "  Written: $means_out"
done

log "  Mean count files in: $MEANS_DIR"

# =============================================================================
# STEP 8 – Fold-change calculations for all group-wise comparisons
# =============================================================================
# Every ordered pair of groups (A vs B) is compared automatically, so
# any new groups added to the sample sheet are handled with no code changes.
#
# fold_change = (mean_A + pseudocount) / (mean_B + pseudocount)
# Pseudocount of 1 avoids division-by-zero and dampens extreme ratios
# when counts are very low. These fold-changes are indicative only.
#
# Output columns (sorted by |fold_change| descending):
#   gene_id  description  mean_A  mean_B  fold_change_A_over_B

log "========== STEP 8: Fold-change calculations =========="

PSEUDOCOUNT=1
n_groups=${#GROUPS[@]}

for (( i=0; i<n_groups; i++ )); do
    for (( j=0; j<n_groups; j++ )); do
        [[ $i -eq $j ]] && continue

        grpA="${GROUPS[$i]}"
        grpB="${GROUPS[$j]}"
        fileA="${MEANS_DIR}/${grpA}_mean_counts.txt"
        fileB="${MEANS_DIR}/${grpB}_mean_counts.txt"

        if [[ ! -f "$fileA" || ! -f "$fileB" ]]; then
            log "  WARNING: Missing means file for ${grpA} or ${grpB}. Skipping."
            continue
        fi

        fc_out="${FC_DIR}/${grpA}_vs_${grpB}_foldchange.txt"
        [[ -f "$fc_out" ]] && { log "  FC file exists: $fc_out – skipping."; continue; }

        log "  Fold-change: ${grpA} vs ${grpB}"

        # AWK two-file join:
        #   File A (NR==FNR): load mean counts and descriptions into arrays
        #   File B          : compute fold-change for each gene
        #   abs_fc column is used only for sorting, then dropped from output
        awk -v pseudo="$PSEUDOCOUNT" '
        BEGIN { OFS="\t" }

        # ---- Reading file A ----
        FNR == NR {
            if ($1 == "gene_id") next
            meanA[$1] = $2+0
            desc[$1]  = (NF>=3) ? $3 : "no_description"
            next
        }

        # ---- Reading file B ----
        {
            if ($1 == "gene_id") next
            gene = $1
            mA   = (gene in meanA) ? meanA[gene] : 0
            mB   = $2+0
            fc   = (mA + pseudo) / (mB + pseudo)
            abs_fc = (fc < 0) ? -fc : fc
            printf "%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n",
                gene,
                (gene in desc) ? desc[gene] : "no_description",
                mA, mB, fc, abs_fc
        }
        ' "$fileA" "$fileB" \
        | sort -t $'\t' -k6,6gr \
        | awk -v a="$grpA" -v b="$grpB" '
            BEGIN {
                OFS="\t"
                print "gene_id","description","mean_"a,"mean_"b,"fold_change_"a"_over_"b
            }
            { print $1,$2,$3,$4,$5 }
        ' > "$fc_out"

        log "  Written: $fc_out"
    done
done

log "  All fold-change files in: $FC_DIR"

# =============================================================================
# PIPELINE COMPLETE
# =============================================================================

log "=========================================="
log "  Pipeline completed successfully!"
log "  Master log : $LOG_FILE"
log "  Outputs    : $OUT_DIR"
log "=========================================="
