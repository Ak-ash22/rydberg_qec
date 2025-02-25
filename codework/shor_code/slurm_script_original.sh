#!/bin/bash

#SBATCH -J q001b
#SBATCH -o logs_auto_runner/%x_%a_log.out
#SBATCH -e logs_auto_runner/%x_%a_log.err
#SBATCH --cpus-per-task=16
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -p epyc-256
#SBATCH --mail-type=FAIL,END
#SBATCH --array=1-100

id=$SLURM_ARRAY_TASK_ID

~/julia-1.11.3/bin/julia main.jl $id
