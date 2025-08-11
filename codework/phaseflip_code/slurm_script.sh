#!/bin/bash

#SBATCH -J s1
#SBATCH -o logs_auto_runner/%x/%A_%a_log.out
#SBATCH -e logs_auto_runner/%x/%A_%a_log.err
#SBATCH --cpus-per-task=4
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -p epyc-768
#SBATCH --mail-type=FAIL,END
#SBATCH --array=0-9

# -------- limit every library to the 4 CPUs you asked for ----------
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK      
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# (optional) avoid depot lock contention
#export JULIA_DEPOT_PATH=$SLURM_TMPDIR/julia_depot

id=$SLURM_ARRAY_TASK_ID

~/julia-1.11.3/bin/julia error_correction_with_storage.jl $id
