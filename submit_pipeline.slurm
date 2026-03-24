#!/usr/bin/env bash
#SBATCH --job-name=TcongoPipeline          # Job name shown in squeue
#SBATCH --output=slurm-%j.out              # Stdout log (%j = job ID)
#SBATCH --error=slurm-%j.err               # Stderr log
#SBATCH --time=08:00:00                    # Max wall time (HH:MM:SS); adjust if needed
#SBATCH --ntasks=1                         # Single task (our pipeline is serial per sample)
#SBATCH --cpus-per-task=4                  # Cores handed to bowtie2/samtools (match THREADS in pipeline)
#SBATCH --mem=16G                          # RAM — bowtie2 indexing and alignment can be memory-hungry
#SBATCH --partition=standard               # Partition/queue name — change to match your HPC (e.g. "compute", "short")
#SBATCH --mail-type=BEGIN,END,FAIL         # Email notifications on job start, finish, or failure
#SBATCH --mail-user=your_email@university.ac.uk   # <-- change to your email

# =============================================================================
# submit_pipeline.slurm
# =============================================================================
# SLURM job submission wrapper for MyFirstPipeline.sh
#
# Usage (from login node):
#   sbatch submit_pipeline.slurm
#
# To monitor the job:
#   squeue -u $USER
#   tail -f slurm-<jobid>.out
#
# To cancel the job:
#   scancel <jobid>
#
# To inspect resource usage after completion:
#   seff <jobid>          (if seff is installed on your cluster)
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,ExitCode
# =============================================================================

set -euo pipefail

# =============================================================================
# SECTION 1 – Print job environment info (helps debugging later)
# =============================================================================

echo "=============================================="
echo "  SLURM Job ID      : $SLURM_JOB_ID"
echo "  Job Name          : $SLURM_JOB_NAME"
echo "  Node              : $SLURMD_NODENAME"
echo "  CPUs allocated    : $SLURM_CPUS_PER_TASK"
echo "  Memory allocated  : ${SLURM_MEM_PER_NODE:-not set}MB"
echo "  Start time        : $(date)"
echo "  Working directory : $(pwd)"
echo "=============================================="

# =============================================================================
# SECTION 2 – Load required modules
# =============================================================================
# On most HPC clusters, tools are not available by default; they are loaded
# via the module system.  Adjust module names to match your cluster's naming.
# Run 'module avail' on your cluster to see what is available.
# If your server has tools installed globally (e.g. the MSc server), you can
# comment out or remove the module load lines.

echo "Loading modules..."

# Uncomment and adjust these as needed for your HPC:
# module purge                        # clear any pre-loaded modules first
# module load fastqc/0.11.9
# module load bowtie2/2.4.4
# module load samtools/1.15
# module load bedtools/2.30.0

# Verify tools are accessible after module loading
echo "Verifying tools in PATH..."
for tool in fastqc bowtie2 bowtie2-build samtools bedtools awk sort; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: '$tool' not found after module load. Check module names."
        exit 1
    fi
    echo "  OK: $tool -> $(command -v $tool)"
done

# =============================================================================
# SECTION 3 – Run the pipeline
# =============================================================================
# The pipeline script is expected to be in the same directory as this SLURM
# script.  We use the directory of this script as a reference so it works
# regardless of where sbatch is called from.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${SCRIPT_DIR}/MyFirstPipeline.sh"

if [[ ! -f "$PIPELINE" ]]; then
    echo "ERROR: Pipeline script not found at: $PIPELINE"
    exit 1
fi

echo "Starting pipeline: $PIPELINE"
echo "----------------------------------------------"

# Run the pipeline; capture exit code without triggering set -e
bash "$PIPELINE"
PIPELINE_EXIT=$?

# =============================================================================
# SECTION 4 – Post-run reporting and failure debugging
# =============================================================================

echo "----------------------------------------------"
echo "  End time: $(date)"

if [[ $PIPELINE_EXIT -ne 0 ]]; then
    echo ""
    echo "  *** PIPELINE FAILED (exit code: $PIPELINE_EXIT) ***"
    echo ""
    echo "  Debugging tips:"
    echo "  1. Check the SLURM stderr log:     cat slurm-${SLURM_JOB_ID}.err"
    echo "  2. Check the pipeline runtime log: ls pipeline_output/pipeline_run_*.log"
    echo "  3. Check per-tool logs:            ls pipeline_output/logs/"
    echo "  4. Check node resource usage:      sacct -j ${SLURM_JOB_ID} --format=JobID,Elapsed,MaxRSS,ExitCode"
    echo "  5. Common causes:"
    echo "       - Out of memory  : increase --mem in this script"
    echo "       - Out of time    : increase --time in this script"
    echo "       - Missing module : run 'module avail' and fix module names in SECTION 2"
    echo "       - Bad file path  : check paths in the CONFIGURATION block of MyFirstPipeline.sh"
    exit $PIPELINE_EXIT
else
    echo "  Pipeline completed successfully!"
    echo ""
    echo "  Resource usage summary:"
    # sacct gives a summary of CPU time, memory, and exit codes for this job
    sacct -j "$SLURM_JOB_ID" \
          --format=JobID,JobName,Elapsed,CPUTime,MaxRSS,ExitCode \
          --noheader 2>/dev/null || echo "  (sacct not available on this system)"
fi

echo "=============================================="
