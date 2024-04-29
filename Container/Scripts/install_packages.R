### Installation Settings ###
# CRAN mirror to get R packages
cran_mirror <- "https://cran.r-project.org/"
proportion_cpu_used <- 1

# List of libraries/packages and version requirements
packages <- list(
  tidyverse = "1.3.1",
  caret = "6.0-94",
  dplyr = "1.1.4",
  magrittr = "2.0.3",
  foreach = "1.5.2",
  doParallel = "1.0.17",
  randomForest = "4.7-1.1",
  kknn = "1.3.1",
  glmnet = "4.1-8",
  yaml = "2.3.8",
  optparse = "1.7.4",
  remotes = "2.5.0"
)

#Install not working above
install.packages("kknn")

cat("Environment Variables Prepared\n")

### Installing Packages ###
# Number of threads/CPU (round up)
ncpus_for_installation <- ceiling(parallel::detectCores() * proportion_cpu_used)

cat("Number of CPU Cores for Installation:", ncpus_for_installation, "\n")

# Install devtools
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", dependencies = TRUE, repos = cran_mirror, Ncpus = ncpus_for_installation)
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Failed to install devtools")
  }
}


# Install specified versions of packages using devtools::install_version
sapply(names(packages), function(pkg_name) {
  tryCatch(
    {
      devtools::install_version(
        pkg_name,
        version = packages[[pkg_name]],
        dependencies = TRUE,
        repos = cran_mirror,
        Ncpus = ncpus_for_installation
      )
    },
    error = function(e) {
      cat("Error installing package '", pkg_name, "': ", conditionMessage(e), "\n")
    }
  )
})

cat("Specific Packages Installed\n")
