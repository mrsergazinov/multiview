#!/bin/bash

Rscript demo2.R sigma1=1 sigma2=1 phi_max=0 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R sigma1=1 sigma2=1 phi_max=0.5 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R sigma1=1 sigma2=1 phi_max=0.8 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R sigma1=2 sigma2=3 phi_max=0 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R sigma1=2 sigma2=3 phi_max=0.5 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R sigma1=2 sigma2=3 phi_max=0.8 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R sigma1=3 sigma2=4 phi_max=0 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R sigma1=3 sigma2=4 phi_max=0.5 no_joint=FALSE no_indiv=FALSE
Rscript demo2.R sigma1=3 sigma2=4 phi_max=0.8 no_joint=FALSE no_indiv=FALSE

Rscript demo2.R sigma1=1 sigma2=1 phi_max=0 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R sigma1=1 sigma2=1 phi_max=0.5 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R sigma1=1 sigma2=1 phi_max=0.8 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R sigma1=2 sigma2=3 phi_max=0 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R sigma1=2 sigma2=3 phi_max=0.5 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R sigma1=2 sigma2=3 phi_max=0.8 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R sigma1=3 sigma2=4 phi_max=0 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R sigma1=3 sigma2=4 phi_max=0.5 no_joint=TRUE no_indiv=FALSE
Rscript demo2.R sigma1=3 sigma2=4 phi_max=0.8 no_joint=TRUE no_indiv=FALSE

Rscript demo2.R sigma1=1 sigma2=1 phi_max=0 no_joint=FALSE no_indiv=TRUE
Rscript demo2.R sigma1=2 sigma2=3 phi_max=0 no_joint=FALSE no_indiv=TRUE
Rscript demo2.R sigma1=3 sigma2=4 phi_max=0 no_joint=FALSE no_indiv=TRUE