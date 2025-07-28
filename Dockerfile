FROM --platform=linux/arm64 rocker/rstudio:4.4.3

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    gdebi-core \
    libcurl4-openssl-dev \
    libssl-dev \
    libclang-dev \
    libxml2-dev \
    libpng-dev \
    libpng-tools \
    psmisc \
    lsb-release \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libssh2-1-dev \
    zlib1g-dev \
    libgmp3-dev \
    libglpk-dev \
    libgsl-dev \
    libmagick++-dev \
    libproj-dev \
    libgdal-dev \
    libudunits2-dev \
    libv8-dev \
    libsqlite3-dev \
    libexpat1-dev \
    default-jdk \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Install renv
RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

# Set working directory
WORKDIR /home/rstudio

# Create user (optional if not using rocker/rstudio)

# Create renv directory with correct permissions
RUN mkdir -p /home/rstudio/renv && \
    chown -R rstudio:rstudio /home/rstudio/renv

# Create .Rprofile to set up renv
RUN echo 'source("renv/activate.R")' > .Rprofile && \
    chown rstudio:rstudio .Rprofile

EXPOSE 8787
# Default command (can change to R or bash depending on your use)
CMD ["/init"]
