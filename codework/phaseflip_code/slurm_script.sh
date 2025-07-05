#!/bin/bash

#SBATCH -J case1
#SBATCH -o logs_auto_runner/%x_%a_log.out
#SBATCH -e logs_auto_runner/%x_%a_log.err
#SBATCH --cpus-per-task=4
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -p epyc-256
#SBATCH --mail-type=FAIL,END
#SBATCH --array=1-1000

# -------- limit every library to the 3 CPUs you asked for ----------
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK      
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# (optional) avoid depot lock contention
#export JULIA_DEPOT_PATH=$SLURM_TMPDIR/julia_depot

id=$SLURM_ARRAY_TASK_ID

~/julia-1.11.3/bin/julia error_detection_and_correction.jl $id
