# Base image for R
ARG BASE_IMAGE="r-base"
ARG BASE_IMAGE_VERSION="4.4.0"
ARG BASE_IMAGE_DIGEST=sha256:2a29b5ba91a84dd60f5706546e6165b292c7e8f4da6c07153c1a0a3a4a01e429
ARG IMAGE_NAME=MnM
ARG IMAGE_VERSION_TAG=1.0.0

# Add platform specification
FROM ${BASE_IMAGE_NAME}@${BASE_IMAGE_DIGEST}

#ARG declarations for image metadata
ARG BASE_IMAGE_NAME
ARG BASE_IMAGE_VERSION
ARG BASE_IMAGE_DIGEST
ARG IMAGE_NAME
ARG IMAGE_VERSION

LABEL org.opencontainers.image.title="${IMAGE_NAME}" \
      org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.authors="PMC Genomicscore" \
      apps.R.source="http://cloud.r-project.org/bin/linux/debian" \
      apps.R.version=${BASE_IMAGE_VERSION} \
      apps.R.licenses="https://spdx.org/licenses/GPL-2.0-or-later"


WORKDIR /opt/
USER root

# System dependencies (apt-get for r-base & libs for specific r package dependencies)
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
      build-essential \
      libssl-dev \
      libcurl4-gnutls-dev \
      libcurl4-openssl-dev \
      libxml2-dev \
      libbz2-dev \
      libpq-dev \
      libsqlite3-dev \
      libcairo2-dev \
      libmariadbd-dev \
      libssh2-1-dev \
      libfontconfig1-dev \
      libfreetype6-dev \
      libpng-dev \
      libtiff5-dev \
      libjpeg-dev \
  && Rscript ./install_packages.R

# Install R packages (CRAN and Bioconductor) using R command block
RUN R -e "install.packages(c('BiocManager', 'dplyr', 'ggplot2'), repos='https://cloud.r-project.org')"

# Install R packages (CRAN packages)
RUN R -e "\
    options(Ncpus = parallel::detectCores(), repos = 'https://cloud.r-project.org'); \
    install.packages(c('BiocManager', 'devtools')); \
    devtools::install_version(c('tidyverse', 'caret', 'dplyr', 'magrittr', 'foreach', 'doParallel', \
        'randomForest', 'kknn', 'glmnet', 'yaml', 'optparse', 'remotes'), \
        version = c('1.3.1', '6.0-94', '1.1.4', '2.0.3', '1.5.2', '1.0.17', '4.7-1.1', \
        '1.3.1', '4.1-8', '2.3.8', '1.7.4', '2.5.0'), dependencies = TRUE); \
    cat('All specified packages installed successfully.'); \
"

# Create and set the working directory
WORKDIR /app

# Create a directory for scripts
#RUN mkdir -p /app/Scripts /app/Inputs /app/MnM /app/Outputs /app/SavedModels

# Copy the "MnM" package from relevant path, not a local directory
#COPY ./MnM /app/MnM


# Install MnM package from local directory (will be changed to install_github when published)
#RUN R -e 'remotes::install_local("/app/MnM/", dependencies = TRUE, force = TRUE)'

# Copy SavedModels data
#COPY ./SavedModels /app/SavedModels

# Copy input data
#COPY ./Inputs /app/Inputs

# Copy all .R files from the local Scripts directory into the container's /app/scripts directory
#COPY ./Scripts/ /app/Scripts/

#Create an entrypoint to run from terminal
#ENTRYPOINT [ "/bin/bash", "-c" ]