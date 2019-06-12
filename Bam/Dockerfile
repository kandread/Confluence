# Use the latest RStudio
FROM rocker/r-ver:3.5.2

# ---------- Core Development ----------

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        apt-utils \
        default-jdk \
        libbz2-dev \
        libcairo2-dev \
        libcurl4-openssl-dev \
        libglu1-mesa-dev \
        libgsl0-dev \
        libhunspell-dev \
        libicu-dev \
        liblzma-dev \
        libmagick++-dev \
        libnlopt-dev \
        libomp-dev \
        libpq-dev \
        libssh2-1-dev \
        libssl-dev \
        libxml2-dev \
        libopenmpi-dev \
        libzmq3-dev \
        mesa-common-dev \
        ssh \
        tk \
        unixodbc-dev \
        libnetcdf-dev \
        && R -e "source('https://bioconductor.org/biocLite.R')"

# ---------- RStan Configuration ----------

# Global site-wide config for building packages
RUN mkdir -p $HOME/.R/ \
    && echo "CXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -flto -ffat-lto-objects  -Wno-unused-local-typedefs \n" >> $HOME/.R/Makevars

# Install ggplot extensions like ggstance and ggrepel
# Install ed, since nloptr needs it to compile
# Install all the dependencies needed by rstanarm and friends
RUN apt-get -y --no-install-recommends install \
    ed \
    clang  \
    ccache \
    && install2.r --error \
        ggstance ggrepel \
        PKI RCurl RJSONIO packrat minqa nloptr matrixStats inline \
        gtools rsconnect shiny \
        xts bayesplot lme4 loo rstantools StanHeaders RcppEigen \
        rstan shinystan rstanarm broom

# ---------- Various R Packages ----------

# Breaking this up in sections due to high chance of packages failing
RUN install2.r --error --deps TRUE devtools formatR selectr caTools remotes DataExplorer
RUN install2.r --error --deps TRUE tidyverse caret GGally outliers hrbrthemes reprex
RUN install2.r --error --deps TRUE lubridate xgboost

RUN R -e "install.packages(c('dplyr','ggplot2','httr','ncdf4'))"
RUN R -e "devtools::install_github('markwh/swotr')"
RUN R -e "devtools::install_github('markwh/bamr')"

# install modules for parallel programming
RUN install2.r snow doSNOW
RUN install2.r foreach iterators
RUN install2.r doParallel doMC doRNG

RUN install2.r future
RUN install2.r future.apply
RUN install2.r doFuture
RUN install2.r future.callr
RUN install2.r furrr

RUN install2.r BatchJobs future.BatchJobs   ## heavy set of dependencies
RUN install2.r batchtools future.batchtools ## heavy set of dependencies
RUN install2.r clustermq


RUN rm -rf /tmp/downloaded_packages/ /tmp/*.rds

COPY . /usr/local/src/app
WORKDIR /usr/local/src/app
# delete the data folder from docker image to save space.
RUN rm -rf /usr/local/src/app/data

# set the unix commands to run the app
CMD ["Rscript","app.R"]
