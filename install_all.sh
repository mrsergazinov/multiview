#!/bin/bash

# Install packages from requirements.txt
python3 -m venv multiview_pylibs
source multiview_pylibs/bin/activate
pip3 install -r install_py_packages.txt

# Install R packages
Rscript install_r_packages.R
