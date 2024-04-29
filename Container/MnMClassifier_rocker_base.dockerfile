# Base image https://hub.docker.com/u/rocker/
ARG BASE_IMAGE="rocker/r-base"
ARG BASE_IMAGE_VERSION="4.3.0"
ARG BASE_IMAGE_DIGEST=sha256:7eb649127377c7d8ff0cce0fb7ac87f0c49642e07d6b48482913befec68381c5
ARG IMAGE_NAME
ARG IMAGE_VERSION_TAG

FROM ${BASE_IMAGE}@${BASE_IMAGE_DIGEST}

USER root

COPY ./Scripts/install_packages.R ./install_packages.R

# System dependencies (apt-get for r-base & libs for specific r package dependencies)
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
      build-essential \
      libssl-dev \
      libcurl4-gnutls-dev \
      libxml2-dev \
      libpq-dev \
      libharfbuzz-dev \
      libfribidi-dev \
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

# Create and set the working directory
WORKDIR /app

# Create a directory for scripts
RUN mkdir -p /app/Scripts /app/Inputs /app/mnm /app/Outputs /app/SavedModels

# Copy the "MnM" package from relevant path, not a local directory
COPY ./MnM /app/MnM

# Via a terminal run, install MnM
#R CMD INSTALL /app/MnM_1.0.0.tar.gz

# Install MnM package from local directory (will be changed to install_github when published)
RUN R -e 'remotes::install_local("/app/MnM/", dependencies = TRUE, force = TRUE)'

# Copy SavedModels data
COPY ./SavedModels /app/SavedModels

# Copy input data
COPY ./Inputs /app/Inputs

# Copy all .R files from the local Scripts directory into the container's /app/scripts directory
COPY ./Scripts/ /app/Scripts/

# Change the execution permission to the wrapper script
RUN chmod +x /app/Scripts/RunningMnM_docker.R

# Run the Wrapper script when the container starts (entry point)
CMD ["Rscript", "/app/Scripts/RunningMnM.R"]

#Create an entrypoint to run from terminal
ENTRYPOINT [ "/bin/bash", "-c" ]