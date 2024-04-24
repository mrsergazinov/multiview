#!/bin/bash

#SBATCH --job-name=jive
#SBATCH --time=2-0
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --partition=long
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --output=out.%j
# module purge
# module load R Python
# bash install_all.sh

Rscript dummy.R