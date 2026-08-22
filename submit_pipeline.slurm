#!/usr/bin/env bash
#SBATCH --job-name=TcongoPipeline
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=standard
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your_email@university.ac.uk

# =============================================================================
# submit_pipeline.slurm
# =============================================================================
# SLURM wrapper that runs the Snakemake workflow inside an Apptainer container.
#
# Usage (from login node):
#   sbatch submit_pipeline.slurm
#
# Monitor:
#   squeue -u $USER
#   tail -f slurm-<jobid>.out
#
# Resource usage after completion:
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,ExitCode
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER="${SCRIPT_DIR}/containers/rnaseq.sif"
DEFINITION="${SCRIPT_DIR}/containers/rnaseq.def"

echo "=============================================="
echo "  Job ID    : $SLURM_JOB_ID"
echo "  Node      : $SLURMD_NODENAME"
echo "  CPUs      : $SLURM_CPUS_PER_TASK"
echo "  Start     : $(date)"
echo "  Directory : $SCRIPT_DIR"
echo "=============================================="

# -- Build Apptainer container if not already built --
if [[ ! -f "$CONTAINER" ]]; then
    echo "Container not found. Building from definition file..."
    apptainer build "$CONTAINER" "$DEFINITION"
    echo "Container built: $CONTAINER"
else
    echo "Container already exists: $CONTAINER"
fi

# -- Run Snakemake inside the Apptainer container --
# --bind mounts the project directory into the container at /pipeline
# so Snakemake can read input files and write outputs

echo "Starting Snakemake workflow..."
echo "----------------------------------------------"

apptainer exec \
    --bind "${SCRIPT_DIR}:/pipeline" \
    "$CONTAINER" \
    snakemake \
        --snakefile "/pipeline/Snakefile" \
        --configfile "/pipeline/config.yaml" \
        --cores "$SLURM_CPUS_PER_TASK" \
        --rerun-incomplete \
        --printshellcmds \
        2>&1

PIPELINE_EXIT=$?

echo "----------------------------------------------"
echo "  End time: $(date)"

if [[ $PIPELINE_EXIT -ne 0 ]]; then
    echo ""
    echo "  *** PIPELINE FAILED (exit code: $PIPELINE_EXIT) ***"
    echo ""
    echo "  Debugging tips:"
    echo "  1. Check SLURM stderr  : cat slurm-${SLURM_JOB_ID}.err"
    echo "  2. Check Snakemake log : cat .snakemake/log/*.log"
    echo "  3. Resource usage      :"
    sacct -j "$SLURM_JOB_ID" \
          --format=JobID,Elapsed,MaxRSS,ExitCode \
          --noheader 2>/dev/null || true
    exit $PIPELINE_EXIT
else
    echo "  Pipeline completed successfully!"
    sacct -j "$SLURM_JOB_ID" \
          --format=JobID,Elapsed,MaxRSS,ExitCode \
          --noheader 2>/dev/null || true
fi

echo "=============================================="
