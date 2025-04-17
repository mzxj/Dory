# Set a default CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

required_packages <- c("optparse", "ggplot2", "reshape2", "dplyr", "data.table", "future", "parallel", "furrr") 

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(new_packages)) {
  install.packages(new_packages, dependencies = TRUE)
}


