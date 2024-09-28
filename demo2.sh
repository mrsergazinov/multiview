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

# Rscript demo_2_views.R snr1=4 snr2=4 phi_max=90 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=4 snr2=4 phi_max=75 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=4 snr2=4 phi_max=60 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=4 snr2=4 phi_max=45 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=75 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=60 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=45 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=1 snr2=1 phi_max=90 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=1 snr2=1 phi_max=75 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=1 snr2=1 phi_max=60 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=1 snr2=1 phi_max=45 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=90 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=75 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=60 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=45 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=90 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=75 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=60 no_joint=FALSE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=45 no_joint=FALSE no_indiv=FALSE

# no joint
# Rscript demo_2_views.R snr1=4 snr2=4 phi_max=90 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=4 snr2=4 phi_max=75 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=4 snr2=4 phi_max=60 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=4 snr2=4 phi_max=45 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=90 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=75 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=60 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=45 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=1 snr2=1 phi_max=90 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=1 snr2=1 phi_max=75 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=1 snr2=1 phi_max=60 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=1 snr2=1 phi_max=45 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=90 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=75 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=60 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=45 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=90 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=75 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=60 no_joint=TRUE no_indiv=FALSE
# Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=45 no_joint=TRUE no_indiv=FALSE

# misspecified ranks
Rscript demo_2_views.R snr1=4 snr2=4 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=4 snr2=4 phi_max=75 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=4 snr2=4 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=4 snr2=4 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=1
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=75 no_joint=FALSE no_indiv=FALSE rank_spec=1
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
# Rscript demo_2_views.R snr1=2 snr2=2 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=1 snr2=1 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=1 snr2=1 phi_max=75 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=1 snr2=1 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=1 snr2=1 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=1
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=75 no_joint=FALSE no_indiv=FALSE rank_spec=1
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
# Rscript demo_2_views.R snr1=0.75 snr2=0.75 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=75 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo_2_views.R snr1=0.5 snr2=0.5 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=1
