# List of packages to install
my_lib_path <- "./multiview_rlibs"
if (!file.exists(my_lib_path)) dir.create(my_lib_path)
.libPaths(my_lib_path)

packages <- c("tidyverse", 
              "r.jive", 
              "RMTstat", 
              "pracma", 
              "Ckmeans.1d.dp", 
              "foreach", 
              "doParallel", 
              "remotes",
              "reticulate")

# Function to check and install packages
install.packages.if.necessary <- function(packages) {
  new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.packages, 
                                             dependencies = TRUE, 
                                             repos = "http://cran.us.r-project.org",
                                             lib = my_lib_path)
}

# Run the installation function
install.packages.if.necessary(packages)

# Install custom packages from Github
remotes::install_github("irinagain/SLIDE")
remotes::install_github("idc9/r_jive")

