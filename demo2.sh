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

Rscript demo2.R snr1=8 snr2=8 phi_max=90 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=60 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=45 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=30 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=90 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=60 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=45 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=30 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=60 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=45 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=30 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=1 snr2=1 phi_max=90 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=1 snr2=1 phi_max=60 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=1 snr2=1 phi_max=45 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=1 snr2=1 phi_max=30 no_joint=FALSE no_indiv=FALSE

# no joint
Rscript demo2.R snr1=8 snr2=8 phi_max=90 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=60 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=45 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=30 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=90 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=60 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=45 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=30 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=90 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=60 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=45 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=30 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=1 snr2=1 phi_max=90 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=1 snr2=1 phi_max=60 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=1 snr2=1 phi_max=45 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=1 snr2=1 phi_max=30 no_joint=TRUE no_indiv=FALSE

# no individuals
Rscript demo2.R snr1=8 snr2=8 phi_max=90 no_joint=FALSE no_indiv=TRUE
Rscript demo2.R snr1=4 snr2=4 phi_max=90 no_joint=FALSE no_indiv=TRUE
Rscript demo2.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=TRUE
Rscript demo2.R snr1=1 snr2=1 phi_max=90 no_joint=FALSE no_indiv=TRUE

# misspecified ranks
Rscript demo2.R snr1=8 snr2=8 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=8 snr2=8 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=8 snr2=8 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=8 snr2=8 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=4 snr2=4 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=4 snr2=4 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=4 snr2=4 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=4 snr2=4 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=2 snr2=2 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=2 snr2=2 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=2 snr2=2 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=1 snr2=1 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=1 snr2=1 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=1 snr2=1 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=1 snr2=1 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=1

Rscript demo2.R snr1=8 snr2=8 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=8 snr2=8 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=8 snr2=8 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=8 snr2=8 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=4 snr2=4 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=4 snr2=4 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=4 snr2=4 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=4 snr2=4 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=2 snr2=2 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=2 snr2=2 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=2 snr2=2 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=1 snr2=1 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=1 snr2=1 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=1 snr2=1 phi_max=45 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=1 snr2=1 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=-1

