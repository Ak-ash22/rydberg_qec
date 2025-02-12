#!/bin/bash

#SBATCH -J q001b
#SBATCH -o logs_auto_runner/%x_log.out
#SBATCH -e logs_auto_runner/%x_log.err
#SBATCH --cpus-per-task=1
#SBATCH --time=200:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH -p epyc-256
#SBATCH --mail-type=FAIL,END
#SBATCH --array=1-1000

id=${SLURM_ARRAY_TASK_ID}

julia main.jl ${id}
