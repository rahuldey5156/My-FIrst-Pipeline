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

set -euo pipefail

echo "=============================================="
echo "  Job ID    : $SLURM_JOB_ID"
echo "  Node      : $SLURMD_NODENAME"
echo "  Start     : $(date)"
echo "=============================================="

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Build the Apptainer container if it doesn't exist yet
if [[ ! -f "${SCRIPT_DIR}/containers/rnaseq.sif" ]]; then
    echo "Building Apptainer container..."
    apptainer build \
        "${SCRIPT_DIR}/containers/rnaseq.sif" \
        "${SCRIPT_DIR}/containers/rnaseq.def"
fi

# Run Snakemake inside the Apptainer container
# --cores        : total cores available to Snakemake
# --use-singularity: tells Snakemake to run rules inside the container
# --singularity-args: passes the container image to use

apptainer exec \
    --bind "${SCRIPT_DIR}" \
    "${SCRIPT_DIR}/containers/rnaseq.sif" \
    snakemake \
        --snakefile "${SCRIPT_DIR}/Snakefile" \
        --configfile "${SCRIPT_DIR}/config.yaml" \
        --cores "$SLURM_CPUS_PER_TASK" \
        --rerun-incomplete \
        --printshellcmds \
        2>&1

EXIT=$?

echo "----------------------------------------------"
echo "  End time: $(date)"

if [[ $EXIT -ne 0 ]]; then
    echo "  *** PIPELINE FAILED (exit code: $EXIT) ***"
    echo "  Check: slurm-${SLURM_JOB_ID}.err"
    sacct -j "$SLURM_JOB_ID" \
          --format=JobID,Elapsed,MaxRSS,ExitCode \
          --noheader 2>/dev/null || true
    exit $EXIT
else
    echo "  Pipeline completed successfully!"
    sacct -j "$SLURM_JOB_ID" \
          --format=JobID,Elapsed,MaxRSS,ExitCode \
          --noheader 2>/dev/null || true
fi
