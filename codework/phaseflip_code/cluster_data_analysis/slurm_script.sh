#!/bin/bash

#SBATCH -J cluster_data_analysis
#SBATCH -o logs_auto_runner/%x_%a_log.out
#SBATCH -e logs_auto_runner/%x_%a_log.err
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -p physik-fleischhauer
#SBATCH --mail-type=FAIL,END

~/julia-1.11.3/bin/julia full_code_analysis.jl
