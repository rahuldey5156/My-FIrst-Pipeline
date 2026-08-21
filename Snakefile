# =============================================================================
# Snakefile
# =============================================================================
# Snakemake workflow for T. congolense RNAi RNAseq analysis.
#
# Rules defined:
#   all             - target rule: defines final outputs Snakemake must produce
#   build_index     - decompress genome and build bowtie2 index
#   fastqc          - run FastQC on each raw fastq file
#   align           - align each sample with bowtie2, convert to sorted BAM
#   count_reads     - bedtools coverage to count reads per gene per sample
#   mean_counts     - per-group mean expression with gene descriptions
#   fold_changes    - all pairwise fold-change comparisons
# =============================================================================

import pandas as pd
import os

# -- Load config --
configfile: "config.yaml"

# -- Parse sample sheet --
# Read the sample sheet and build sample metadata.
# Expected columns: SampleName SampleType Replicate Time Treatment End1 End2
samples_df = pd.read_csv(
    config["sample_sheet"],
    sep="\t",
    comment="#"
)

# Build a dictionary: sample_name -> {r1, r2, group}
SAMPLES = {}
for _, row in samples_df.iterrows():
    name  = row["SampleName"]
    group = f"{row['SampleType']}_{row['Treatment']}_T{row['Time']}"
    SAMPLES[name] = {
        "r1":    os.path.join(config["fastq_dir"], row["End1"]),
        "r2":    os.path.join(config["fastq_dir"], row["End2"]),
        "group": group
    }

# Get all unique groups
GROUPS = list(set(s["group"] for s in SAMPLES.values()))

# Get all pairwise group combinations for fold-change
FC_PAIRS = [
    f"{a}_vs_{b}"
    for a in GROUPS
    for b in GROUPS
    if a != b
]

# =============================================================================
# RULE: all
# =============================================================================
# The 'all' rule defines the final outputs the pipeline must produce.
# Snakemake works backwards from these targets to figure out which
# rules to run and in what order.

rule all:
    input:
        # FastQC reports for every read file
        expand(
            "{out}/fastqc/{sample}_{end}_fastqc.zip",
            out=config["out_dir"],
            sample=list(SAMPLES.keys()),
            end=["1", "2"]
        ),
        # Sorted BAM index for every sample
        expand(
            "{out}/bam/{sample}.sorted.bam.bai",
            out=config["out_dir"],
            sample=list(SAMPLES.keys())
        ),
        # Count files for every sample
        expand(
            "{out}/counts/{sample}_counts.txt",
            out=config["out_dir"],
            sample=list(SAMPLES.keys())
        ),
        # Mean count files for every group
        expand(
            "{out}/means/{group}_mean_counts.txt",
            out=config["out_dir"],
            group=GROUPS
        ),
        # Fold-change files for every pair
        expand(
            "{out}/fold_changes/{pair}_foldchange.txt",
            out=config["out_dir"],
            pair=FC_PAIRS
        )

# =============================================================================
# RULE: build_index
# =============================================================================
# Decompress the gzipped genome FASTA and build a bowtie2 index.
# Only runs once — Snakemake skips it if the index already exists.

rule build_index:
    input:
        genome = config["genome_fasta"]
    output:
        # bowtie2-build creates 6 index files; we track the .1.bt2 as sentinel
        index = config["index_prefix"] + ".1.bt2",
        genome_cat = config["index_prefix"] + "_genome.fa"
    threads:
        config["threads"]
    log:
        config["out_dir"] + "/logs/bowtie2_build.log"
    shell:
        """
        # Decompress genome into a single FASTA
        gunzip -c {input.genome} > {output.genome_cat}

        # Build bowtie2 index
        bowtie2-build \
            --threads {threads} \
            {output.genome_cat} \
            {config[index_prefix]} \
            > {log} 2>&1
        """

# =============================================================================
# RULE: fastqc
# =============================================================================
# Run FastQC on each individual fastq.gz file.
# Snakemake runs this rule independently for each sample and each end,
# which means it can run them in parallel if resources allow.

rule fastqc:
    input:
        r1 = lambda wc: SAMPLES[wc.sample]["r1"],
        r2 = lambda wc: SAMPLES[wc.sample]["r2"]
    output:
        zip_r1 = config["out_dir"] + "/fastqc/{sample}_1_fastqc.zip",
        zip_r2 = config["out_dir"] + "/fastqc/{sample}_2_fastqc.zip"
    threads: 2
    log:
        config["out_dir"] + "/logs/{sample}_fastqc.log"
    shell:
        """
        fastqc \
            --outdir {config[out_dir]}/fastqc \
            --threads {threads} \
            --quiet \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """

# =============================================================================
# RULE: align
# =============================================================================
# Align paired-end reads with bowtie2.
# SAM output is piped directly to samtools view (SAM->BAM), then
# samtools sort, avoiding any temporary files on disk.
# samtools index creates the .bai file needed by bedtools downstream.

rule align:
    input:
        r1    = lambda wc: SAMPLES[wc.sample]["r1"],
        r2    = lambda wc: SAMPLES[wc.sample]["r2"],
        index = config["index_prefix"] + ".1.bt2"
    output:
        bam = config["out_dir"] + "/bam/{sample}.sorted.bam",
        bai = config["out_dir"] + "/bam/{sample}.sorted.bam.bai"
    threads:
        config["threads"]
    log:
        config["out_dir"] + "/logs/{sample}_bowtie2.log"
    shell:
        """
        bowtie2 \
            -x {config[index_prefix]} \
            -1 {input.r1} \
            -2 {input.r2} \
            --no-unal \
            --no-mixed \
            --no-discordant \
            -p {threads} \
            2> {log} \
        | samtools view -bS -F 4 - \
        | samtools sort -@ {threads} -o {output.bam} -

        samtools index {output.bam}
        """

# =============================================================================
# RULE: count_reads
# =============================================================================
# Use bedtools coverage to count reads overlapping each gene in the BED file.
# Extracts gene_id (col 4) and count (last col) into a clean two-column file.

rule count_reads:
    input:
        bam = config["out_dir"] + "/bam/{sample}.sorted.bam",
        bai = config["out_dir"] + "/bam/{sample}.sorted.bam.bai",
        bed = config["bed_file"]
    output:
        counts = config["out_dir"] + "/counts/{sample}_counts.txt"
    log:
        config["out_dir"] + "/logs/{sample}_bedtools.log"
    shell:
        """
        bedtools coverage \
            -counts \
            -a {input.bed} \
            -b {input.bam} \
            2> {log} \
        | awk 'BEGIN{{OFS="\t"}} {{print $4, $NF}}' \
        > {output.counts}
        """

# =============================================================================
# RULE: mean_counts
# =============================================================================
# For each group, collect all replicate count files and compute per-gene means.
# Calls an external Python helper script for the calculation and description join.

rule mean_counts:
    input:
        # Collect count files for all samples belonging to this group
        counts = lambda wc: [
            config["out_dir"] + f"/counts/{s}_counts.txt"
            for s, info in SAMPLES.items()
            if info["group"] == wc.group
        ],
        bed = config["bed_file"]
    output:
        means = config["out_dir"] + "/means/{group}_mean_counts.txt"
    log:
        config["out_dir"] + "/logs/{group}_means.log"
    script:
        "scripts/mean_counts.py"

# =============================================================================
# RULE: fold_changes
# =============================================================================
# Compute fold-changes for one group-pair (encoded in the wildcard as A_vs_B).
# Calls an external Python helper script.

rule fold_changes:
    input:
        # Parse the pair name to find the two means files needed
        fileA = lambda wc: config["out_dir"] + "/means/" + wc.pair.split("_vs_")[0] + "_mean_counts.txt",
        fileB = lambda wc: config["out_dir"] + "/means/" + wc.pair.split("_vs_")[1] + "_mean_counts.txt"
    output:
        fc = config["out_dir"] + "/fold_changes/{pair}_foldchange.txt"
    params:
        pseudocount = config["pseudocount"]
    log:
        config["out_dir"] + "/logs/{pair}_fc.log"
    script:
        "scripts/fold_changes.py"
