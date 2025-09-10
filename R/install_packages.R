# Set a default CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

required_packages <- c("optparse", "ggplot2", "reshape2", "dplyr", "data.table", "future", "parallel", "furrr") 
# Find which are not yet installed
installed <- installed.packages()[, "Package"]
new_packages <- required_packages[!(required_packages %in% installed)]
existing_packages <- required_packages[required_packages %in% installed]

# Print summary
message(length(required_packages), " packages are required: ",
        paste(required_packages, collapse = ", "))

if (length(existing_packages) > 0) {
  message(length(existing_packages), " packages already exist: ",
          paste(existing_packages, collapse = ", "))
}

if (length(new_packages) > 0) {
  message(length(new_packages), " packages will now be installed: ",
          paste(new_packages, collapse = ", "))
  install.packages(new_packages, dependencies = TRUE)
}

# Re-check after installation
final_installed <- installed.packages()[, "Package"]
missing_after <- required_packages[!(required_packages %in% final_installed)]

if (length(missing_after) == 0) {
  message("All required packages are now installed successfully.")
} else {
  message("Some packages are still missing: ",
          paste(missing_after, collapse = ", "))
}
