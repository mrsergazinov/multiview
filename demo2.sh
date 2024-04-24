#!/bin/bash

Rscript demo2.R snr1=8 snr2=8 phi_max=90 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=60 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=30 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=90 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=60 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=30 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=60 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=30 no_joint=FALSE no_indiv=FALSE

# no joint
Rscript demo2.R snr1=8 snr2=8 phi_max=90 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=60 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=30 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=90 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=60 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=30 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=90 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=60 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=30 no_joint=TRUE no_indiv=FALSE

# no individuals
Rscript demo2.R snr1=8 snr2=8 phi_max=90 no_joint=FALSE no_indiv=TRUE
Rscript demo2.R snr1=4 snr2=4 phi_max=90 no_joint=FALSE no_indiv=TRUE
Rscript demo2.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=TRUE

# misspecified ranks
Rscript demo2.R snr1=8 snr2=8 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=8 snr2=8 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=8 snr2=8 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=4 snr2=4 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=4 snr2=4 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=4 snr2=4 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=2 snr2=2 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=1
Rscript demo2.R snr1=2 snr2=2 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=1

Rscript demo2.R snr1=8 snr2=8 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=8 snr2=8 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=8 snr2=8 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=4 snr2=4 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=4 snr2=4 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=4 snr2=4 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=2 snr2=2 phi_max=90 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=2 snr2=2 phi_max=60 no_joint=FALSE no_indiv=FALSE rank_spec=-1
Rscript demo2.R snr1=2 snr2=2 phi_max=30 no_joint=FALSE no_indiv=FALSE rank_spec=-1

