#!/bin/bash

Rscript demo2.R snr1=8 snr2=8 phi_max=0 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=0.5 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=0.8 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=0 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=0.5 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=0.8 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=0 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=0.5 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=0.8 no_joint=FALSE no_indiv=FALSE

Rscript demo2.R snr1=8 snr2=8 phi_max=0 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=0.5 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=8 snr2=8 phi_max=0.8 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=0 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=0.5 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=4 snr2=4 phi_max=0.8 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=0 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=0.5 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R snr1=2 snr2=2 phi_max=0.8 no_joint=TRUE no_indiv=FALSE

Rscript demo2.R snr1=8 snr2=8 phi_max=0 no_joint=FALSE no_indiv=TRUE
Rscript demo2.R snr1=4 snr2=4 phi_max=0 no_joint=FALSE no_indiv=TRUE
Rscript demo2.R snr1=2 snr2=2 phi_max=0 no_joint=FALSE no_indiv=TRUE