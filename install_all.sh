#!/bin/bash

# Install packages from requirements.txt
pip install -r install_py_packages.txt

# Install R packages
Rscript install_r_packages.R
