# MyFirstPipeline

A reproducible RNAseq analysis pipeline for *Trypanosoma congolense* RNAi
knock-down experiments. Available in two versions — a standalone bash/AWK
script and a Snakemake workflow with Apptainer containerisation.

## Repository

https://github.com/rahuldey5156/My-FIrst-Pipeline

---

## Two versions of this pipeline

### Version 1 — `MyFirstPipeline.sh`
Original pipeline written in pure bash and AWK. Runs as a single
script with no external workflow manager required. Written as part
of an MSc Bioinformatics coursework assignment.

```bash
bash MyFirstPipeline.sh
```

### Version 2 — `Snakefile` (extended version)
The same pipeline rebuilt as a Snakemake workflow, containerised with
Apptainer and deployable on HPC clusters via SLURM. All pipeline logic
remains in pure bash and AWK inside Snakemake `shell:` blocks.
Snakemake manages rule dependencies, execution order, and parallelisation.

```bash
# Local run
snakemake --configfile config.yaml --cores 4

# HPC cluster
sbatch submit_pipeline.slurm
```

---

## Pipeline steps

| Step | Description |
|------|-------------|
| 1 | Sanity checks — tools, paths, sample sheet, fastq files |
| 2 | FastQC on all raw paired-end fastq.gz files |
| 3 | Parse FastQC summaries → PASS/WARN/FAIL report |
| 4 | Build bowtie2 index from *T. congolense* genome (once only) |
| 5 | Align reads with bowtie2 → sorted indexed BAM |
| 6 | Count reads per gene with bedtools coverage |
| 7 | Per-group mean counts with gene descriptions |
| 8 | All pairwise fold-change comparisons sorted by absolute FC |

---

## Repository structure

```
My-FIrst-Pipeline/
├── MyFirstPipeline.sh        Version 1: standalone bash/AWK pipeline
├── Snakefile                 Version 2: Snakemake workflow
├── config.yaml               Configuration for Snakemake version
├── submit_pipeline.slurm     SLURM job submission wrapper
├── containers/
│   └── rnaseq.def            Apptainer container definition
└── .gitignore
```

---

## Sample sheet format

`Tco2.fqfiles` — tab-delimited, 7 columns:

| SampleName | SampleType | Replicate | Time | Treatment | End1 | End2 |
|------------|------------|-----------|------|-----------|------|------|
| Tco-106 | WT | 1 | 0 | uninduced | Tco-106_1.fq.gz | Tco-106_2.fq.gz |

Group labels are built automatically as `SampleType_Treatment_TTime`
e.g. `Clone1_induced_T24`, `WT_uninduced_T0`

---

## Requirements

### Version 1 (bash script)
The following tools must be available in PATH:

- `fastqc`
- `bowtie2` and `bowtie2-build`
- `samtools`
- `bedtools`
- `awk`, `sort` (standard Unix)

### Version 2 (Snakemake)
Only two things needed on the host:

- `snakemake`
- `apptainer`

All bioinformatics tools (FastQC, bowtie2, samtools, bedtools) are
installed inside the Apptainer container — nothing else required.

---

## Setup for Version 1 (bash script)

**Step 1 — Clone the repository**
```bash
git clone https://github.com/rahuldey5156/My-FIrst-Pipeline.git
cd My-FIrst-Pipeline
```

**Step 2 — Edit paths in the CONFIGURATION block**
```bash
nano MyFirstPipeline.sh
# Update BASE_DIR, FASTQ_DIR, SAMPLE_SHEET, GENOME_DIR, BED_FILE
```

**Step 3 — Run**
```bash
bash MyFirstPipeline.sh
```

---

## Setup for Version 2 (Snakemake)

**Step 1 — Clone the repository**
```bash
git clone https://github.com/rahuldey5156/My-FIrst-Pipeline.git
cd My-FIrst-Pipeline
```

**Step 2 — Edit paths in config.yaml**
```bash
nano config.yaml
```

**Step 3 — Build the Apptainer container (once only)**
```bash
apptainer build containers/rnaseq.sif containers/rnaseq.def
```

**Step 4 — Dry run to check everything**
```bash
snakemake --configfile config.yaml --cores 4 --dry-run
```

**Step 5 — Full run**
```bash
# Local
snakemake --configfile config.yaml --cores 4

# HPC cluster with SLURM
sbatch submit_pipeline.slurm
```

---

## Monitoring a SLURM job

```bash
# Check job status
squeue -u $USER

# Watch live output
tail -f slurm-<jobid>.out

# Check resource usage after completion
sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,ExitCode
```

---

## Output structure

```
pipeline_output/
├── fastqc/           FastQC HTML and zip reports per sample
├── QC_summary.txt    Combined PASS/WARN/FAIL table across all samples
├── bowtie2_index/    Genome index files
├── bam/              Sorted indexed BAM files per sample
├── counts/           Per-sample read counts per gene
├── means/            Per-group mean counts with gene descriptions
├── fold_changes/     All pairwise fold-change tables sorted by |FC|
└── logs/             Per-tool log files
```

---

## Experimental design

The pipeline was developed for a *Trypanosoma congolense* RNAi knock-down
experiment in the IL3000 laboratory strain. The experimental design includes:

- **Time points:** T=0h, T=24h, T=48h
- **Sample types:** Wild Type (WT), Clone1, Clone2 (two independent RNAi lines)
- **Conditions:** induced (tetracycline-treated) and uninduced
- **Replicates:** 3 or more biological replicates per group

Group labels are constructed as `SampleType_Treatment_TTime`, giving
comparisons such as `Clone1_induced_T24_vs_WT_uninduced_T0`.

---

## QC note

All samples show a FastQC FAIL for *Per base sequence content*. This is
completely expected in RNAseq data due to random hexamer priming during
library preparation, which causes non-random base composition at the start
of reads. It does not indicate a data quality problem and all samples
proceed through the pipeline normally.

---

## Debugging failures

**Version 1 (bash script)**
```bash
# Check the master log
cat pipeline_output/pipeline_run_*.log

# Check individual tool logs
ls pipeline_output/logs/
cat pipeline_output/logs/bowtie2_build.log
cat pipeline_output/logs/<sample>_bowtie2.log
```

**Version 2 (Snakemake)**
```bash
# Check Snakemake log
cat .snakemake/log/*.log

# Check individual tool logs
ls pipeline_output/logs/

# Check SLURM resource usage
sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,ExitCode
```

**Common causes of failure**

| Problem | Fix |
|---------|-----|
| Out of memory | Increase `--mem` in `submit_pipeline.slurm` |
| Out of time | Increase `--time` in `submit_pipeline.slurm` |
| Tool not found | Check `module avail` or rebuild Apptainer container |
| Wrong paths | Edit `config.yaml` or the CONFIGURATION block in `MyFirstPipeline.sh` |
| Sample sheet format | Ensure 7 tab-separated columns with no trailing spaces |

---

## Author

Rahul Dey
MSc Bioinformatics, University of Edinburgh
