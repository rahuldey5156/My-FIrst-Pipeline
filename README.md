# MyFirstPipeline

An RNAseq analysis pipeline for *Trypanosoma congolense* RNAi knock-down experiments,
written in bash/awk for the BPSM MSc course.

## Repository

https://github.com/rahuldey5156/My-FIrst-Pipeline

## Pipeline steps

| Step | Description |
|------|-------------|
| 1 | Sanity checks — tools, files, sample sheet validation |
| 2 | FastQC on all raw paired-end fastq.gz files |
| 3 | Parse FastQC summaries — PASS/WARN/FAIL report |
| 4 | Build bowtie2 index from T. congolense genome |
| 5 | Align reads with bowtie2 → sorted indexed BAM (samtools) |
| 6 | Count reads per gene with bedtools coverage |
| 7 | Per-group mean counts with gene descriptions |
| 8 | Pairwise fold-change calculations (all group combinations) |

## Requirements

The following tools must be available in `PATH` (or loaded via modules):

- `fastqc`
- `bowtie2` + `bowtie2-build`
- `samtools`
- `bedtools`
- `awk`, `sort` (standard Unix)

## Running on the MSc server (directly)
```bash
chmod +x MyFirstPipeline.sh
bash MyFirstPipeline.sh
```

## Running on an HPC cluster with SLURM
```bash
# Edit the email address and partition name in submit_pipeline.slurm first, then:
sbatch submit_pipeline.slurm

# Monitor the job
squeue -u $USER

# Stream the live log
tail -f slurm-<jobid>.out

# Check resource usage after completion
sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,ExitCode
```

## Input files expected

| File | Description |
|------|-------------|
| `Tco2.fqfiles` | Tab-delimited sample sheet: sample, R1, R2, group |
| `fastq/` | Directory of paired-end `*.fq.gz` files |
| `Tcongo_genome/` | Directory of genome FASTA files (`.fa` / `.fasta`) |
| `TriTrypDB-46_TcongolenseIL3000_2019.bed` | Gene annotation BED file |

## Output structure
```
pipeline_output/
├── fastqc/               FastQC HTML and zip reports
├── QC_summary.txt        Combined PASS/WARN/FAIL summary
├── bowtie2_index/        Genome index files
├── bam/                  Sorted, indexed BAM files
├── counts/               Per-sample read counts per gene
├── means/                Per-group mean counts with descriptions
├── fold_changes/         All pairwise fold-change tables
└── logs/                 Per-tool log files + master pipeline log
```

## Debugging failures

If the pipeline fails:
1. Check `slurm-<jobid>.err` (SLURM runs) or terminal output (direct runs)
2. Check `pipeline_output/pipeline_run_*.log` for the step that failed
3. Check `pipeline_output/logs/` for the specific tool that failed
4. Common issues:
   - **Out of memory** — increase `--mem` in `submit_pipeline.slurm`
   - **Out of time** — increase `--time` in `submit_pipeline.slurm`
   - **Missing tool** — check `module avail` and update module names
   - **Wrong paths** — edit the `CONFIGURATION` block in `MyFirstPipeline.sh`
