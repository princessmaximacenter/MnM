# Base image for R
ARG BASE_IMAGE_NAME="r-base"
ARG BASE_IMAGE_VERSION=4.4.0
ARG BASE_IMAGE_DIGEST=sha256:2a29b5ba91a84dd60f5706546e6165b292c7e8f4da6c07153c1a0a3a4a01e429
ARG IMAGE_NAME=MnM
ARG IMAGE_VERSION_TAG=1.0.0

FROM ${BASE_IMAGE_NAME}@${BASE_IMAGE_DIGEST}

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

# Install required system dependencies
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y \
      build-essential \
      libssl-dev \
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
      libgit2-dev \
      libharfbuzz-dev \
      libfribidi-dev \
      pkg-config && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e " \
    options(Ncpus = parallel::detectCores(), repos = 'https://cloud.r-project.org'); \
    install.packages(c('BiocManager', 'devtools'), dependencies = TRUE); \
    pkgs <- c('tidyverse', 'caret', 'dplyr', 'magrittr', 'foreach', 'doParallel', \
              'randomForest', 'kknn', 'glmnet', 'yaml', 'optparse', 'remotes'); \
    versions <- c('1.3.1', '6.0-94', '1.1.4', '2.0.3', '1.5.2', '1.0.17', \
                  '4.7-1.1', '1.3.1', '4.1-8', '2.3.8', '1.7.4', '2.5.0'); \
    for (i in seq_along(pkgs)) { \
        tryCatch({ \
            devtools::install_version(pkgs[i], version = versions[i], dependencies = TRUE); \
        }, error = function(e) { \
            cat(sprintf('Failed to install %s. Error: %s\n', pkgs[i], e$message)); \
        }); \
    }; \
    cat('All specified packages installed successfully.\n'); \
"

# Set the working directory
WORKDIR /app

# Install MnM package from github (https://github.com/princessmaximacenter/MnM), specifically dev branch 
RUN R -e " \
    library(remotes); \
    remotes::install_github('princessmaximacenter/MnM', ref = 'dev', dependencies = TRUE, force = TRUE); \
"

# Create a directory for scripts
#RUN mkdir -p /app/Scripts /app/Inputs /app/MnM /app/Outputs /app/SavedModels

#Create an entrypoint to run from terminal
ENTRYPOINT [ "/bin/bash"]

# Copy all .R files from the local Scripts directory into the container's /app/scripts directory
#COPY ./Scripts/ /app/Scripts/