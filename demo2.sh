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

Rscript demo_2_views.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=2 snr2=2 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo_2_views.R snr1=2 snr2=2 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo_2_views.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=0
Rscript demo_2_views.R snr1=2 snr2=2 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=0
Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=0
Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=0