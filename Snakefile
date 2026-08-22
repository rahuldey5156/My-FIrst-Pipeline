# =============================================================================
# Snakefile
# =============================================================================
# Snakemake workflow for T. congolense RNAi RNAseq analysis.
# All pipeline logic is written in pure bash/AWK inside shell: blocks.
# Python is only used at the top to parse the sample sheet and build
# the sample/group lists for Snakemake to work with.
#
# Usage:
#   snakemake --configfile config.yaml --cores 4 --dry-run
#   snakemake --configfile config.yaml --cores 4
#   sbatch submit_pipeline.slurm
# =============================================================================

import csv
import os

configfile: "config.yaml"

# =============================================================================
# PARSE SAMPLE SHEET
# =============================================================================

SAMPLES = {}
with open(config["sample_sheet"]) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        if row["SampleName"].startswith("#"):
            continue
        if row["SampleName"] == "SampleName":
            continue
        name  = row["SampleName"]
        group = "{}_{}_{}" .format(
            row["SampleType"],
            row["Treatment"],
            "T" + str(row["Time"])
        )
        SAMPLES[name] = {
            "r1":    row["End1"],
            "r2":    row["End2"],
            "group": group
        }

GROUPS = sorted(set(s["group"] for s in SAMPLES.values()))

FC_PAIRS = [
    "{}_vs_{}".format(a, b)
    for a in GROUPS
    for b in GROUPS
    if a != b
]

# =============================================================================
# RULE: all
# =============================================================================

rule all:
    input:
        expand(
            os.path.join(config["out_dir"], "fastqc", "{sample}_{end}_fastqc.zip"),
            sample=list(SAMPLES.keys()),
            end=["1", "2"]
        ),
        expand(
            os.path.join(config["out_dir"], "bam", "{sample}.sorted.bam.bai"),
            sample=list(SAMPLES.keys())
        ),
        expand(
            os.path.join(config["out_dir"], "counts", "{sample}_counts.txt"),
            sample=list(SAMPLES.keys())
        ),
        expand(
            os.path.join(config["out_dir"], "means", "{group}_mean_counts.txt"),
            group=GROUPS
        ),
        expand(
            os.path.join(config["out_dir"], "fold_changes", "{pair}_foldchange.txt"),
            pair=FC_PAIRS
        ),
        os.path.join(config["out_dir"], "QC_summary.txt")

# =============================================================================
# RULE: build_index
# =============================================================================

rule build_index:
    input:
        genome = config["genome_fasta"]
    output:
        sentinel   = config["index_prefix"] + ".1.bt2",
        genome_cat = config["index_prefix"] + "_genome.fa"
    threads:
        config["threads"]
    log:
        os.path.join(config["out_dir"], "logs", "bowtie2_build.log")
    shell:
        """
        gunzip -c {input.genome} > {output.genome_cat}

        bowtie2-build \
            --threads {threads} \
            {output.genome_cat} \
            {config[index_prefix]} \
            > {log} 2>&1
        """

# =============================================================================
# RULE: fastqc
# =============================================================================

rule fastqc:
    input:
        r1 = lambda wc: os.path.join(
            config["fastq_dir"], SAMPLES[wc.sample]["r1"]),
        r2 = lambda wc: os.path.join(
            config["fastq_dir"], SAMPLES[wc.sample]["r2"])
    output:
        zip_r1 = os.path.join(
            config["out_dir"], "fastqc", "{sample}_1_fastqc.zip"),
        zip_r2 = os.path.join(
            config["out_dir"], "fastqc", "{sample}_2_fastqc.zip")
    threads: 2
    log:
        os.path.join(config["out_dir"], "logs", "{sample}_fastqc.log")
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
# RULE: qc_summary
# =============================================================================

rule qc_summary:
    input:
        expand(
            os.path.join(config["out_dir"], "fastqc", "{sample}_{end}_fastqc.zip"),
            sample=list(SAMPLES.keys()),
            end=["1", "2"]
        )
    output:
        summary = os.path.join(config["out_dir"], "QC_summary.txt")
    shell:
        """
        printf "%-50s %-40s %s\n" "Sample" "Module" "Result" \
            > {output.summary}
        printf "%100s\n" | tr ' ' '=' >> {output.summary}

        for zip in {config[out_dir]}/fastqc/*_fastqc.zip; do
            [ -f "$zip" ] || continue
            tmpdir=$(mktemp -d)
            unzip -q "$zip" "*/summary.txt" -d "$tmpdir" 2>/dev/null || {{
                echo "WARNING: Could not extract summary from $zip"
                rm -rf "$tmpdir"
                continue
            }}
            summary_file=$(find "$tmpdir" -name "summary.txt" | head -1)
            while IFS=$'\t' read -r result module filename; do
                printf "%-50s %-40s %s\n" \
                    "$filename" "$module" "$result" >> {output.summary}
            done < "$summary_file"
            rm -rf "$tmpdir"
        done
        """

# =============================================================================
# RULE: align
# =============================================================================

rule align:
    input:
        r1       = lambda wc: os.path.join(
            config["fastq_dir"], SAMPLES[wc.sample]["r1"]),
        r2       = lambda wc: os.path.join(
            config["fastq_dir"], SAMPLES[wc.sample]["r2"]),
        sentinel = config["index_prefix"] + ".1.bt2"
    output:
        bam = os.path.join(
            config["out_dir"], "bam", "{sample}.sorted.bam"),
        bai = os.path.join(
            config["out_dir"], "bam", "{sample}.sorted.bam.bai")
    threads:
        config["threads"]
    log:
        os.path.join(config["out_dir"], "logs", "{sample}_bowtie2.log")
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

rule count_reads:
    input:
        bam = os.path.join(
            config["out_dir"], "bam", "{sample}.sorted.bam"),
        bai = os.path.join(
            config["out_dir"], "bam", "{sample}.sorted.bam.bai"),
        bed = config["bed_file"]
    output:
        counts = os.path.join(
            config["out_dir"], "counts", "{sample}_counts.txt")
    log:
        os.path.join(config["out_dir"], "logs", "{sample}_bedtools.log")
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

rule mean_counts:
    input:
        counts = lambda wc: [
            os.path.join(config["out_dir"], "counts", "{}_counts.txt".format(s))
            for s, info in SAMPLES.items()
            if info["group"] == wc.group
        ],
        bed = config["bed_file"]
    output:
        means = os.path.join(
            config["out_dir"], "means", "{group}_mean_counts.txt")
    shell:
        """
        gene_desc=$(mktemp)
        awk 'BEGIN{{OFS="\t"}} NF>=5 {{print $4, $5}}' {input.bed} \
            | sort -u > "$gene_desc"

        means_tmp=$(mktemp)
        awk '
        BEGIN {{ OFS="\t" }}
        {{ sum[$1] += $2+0 }}
        END {{
            n = ARGC - 1
            for (gene in sum)
                printf "%s\t%.4f\n", gene, sum[gene]/n
        }}
        ' {input.counts} | sort -k1,1 > "$means_tmp"

        join -t $'\t' -1 1 -2 1 -a 1 "$means_tmp" "$gene_desc" \
        | awk '
            BEGIN {{ OFS="\t"; print "gene_id","mean_count","description" }}
            {{ print $1, $2, (NF>=3 ? $3 : "no_description") }}
        ' > {output.means}

        rm -f "$means_tmp" "$gene_desc"
        """

# =============================================================================
# RULE: fold_changes
# =============================================================================

rule fold_changes:
    input:
        fileA = lambda wc: os.path.join(
            config["out_dir"], "means",
            wc.pair.split("_vs_")[0] + "_mean_counts.txt"),
        fileB = lambda wc: os.path.join(
            config["out_dir"], "means",
            wc.pair.split("_vs_")[1] + "_mean_counts.txt")
    output:
        fc = os.path.join(
            config["out_dir"], "fold_changes", "{pair}_foldchange.txt")
    params:
        pseudocount = config["pseudocount"]
    shell:
        """
        pair="{wildcards.pair}"
        grpA="${{pair%%_vs_*}}"
        grpB="${{pair##*_vs_}}"

        awk -v pseudo="{params.pseudocount}" \
            -v grpA="$grpA" \
            -v grpB="$grpB" '
        BEGIN {{ OFS="\t" }}

        FNR == NR {{
            if ($1 == "gene_id") next
            meanA[$1] = $2+0
            desc[$1]  = (NF>=3) ? $3 : "no_description"
            next
        }}

        {{
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
        }}
        ' {input.fileA} {input.fileB} \
        | sort -t $'\t' -k6,6gr \
        | awk -v a="$grpA" -v b="$grpB" '
            BEGIN {{
                OFS="\t"
                print "gene_id","description","mean_"a,"mean_"b,"fold_change_"a"_over_"b
            }}
            {{ print $1,$2,$3,$4,$5 }}
        ' > {output.fc}
        """
