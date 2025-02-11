#!/bin/bash

#SBATCH -J q001b
#SBATCH -o logs_auto_runner/%x_log.out
#SBATCH -e logs_auto_runner/%x_log.err
#SBATCH --cpus-per-task=1
#SBATCH --time=200:00:00
#SBATCH --mem=4G
#SBATCH -p epyc-256
#SBATCH --mail-type=FAIL,END
#SBATCH --array=0-1

id=${SLURM_ARRAY_TASK_ID}

module load julia
julia main.jl ${id}
