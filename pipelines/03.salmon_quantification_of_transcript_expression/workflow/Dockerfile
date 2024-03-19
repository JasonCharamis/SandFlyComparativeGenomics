FROM condaforge/mambaforge:latest
LABEL maintainer = "Jason Charamis" 
LABEL contact    = "jason.charamis@gmail.com"
LABEL build_date = "2023-11-30"
LABEL version    = "v.0.0.1-dev"

## Dockerfile to build Docker image for this Snakemake workflow.
FROM ubuntu:20.04

# Set environment variables for non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe

# Step 1: Install essential dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        git \
        python3-pip \
	wget

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Set conda on the PATH
ENV PATH /opt/conda/bin:$PATH
RUN conda init --all

# Step 2: Install R without prompts
RUN apt-get install -y --no-install-recommends \
        r-base \
        r-base-dev \
        && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

ENV PATH /usr/lib/R/bin:$PATH

# Step 2: Retrieve conda environments

# Conda environment:
#   source: envs/salmon.yaml
#   name: salmon
#   channels:
#       - bioconda
#   dependencies:
#       - salmon

RUN conda install -c bioconda salmon &&\
    conda install -c conda-forge r-essentials

# Install BiocManager
RUN R -e 'install.packages("BiocManager", repos="https://cran.rstudio.com")'

# Install required R packages using BiocManager
RUN R -e 'BiocManager::install(c("limma", "voom", "edgeR"))'

# Install other required R packages from CRAN
RUN R -e 'install.packages(c("dplyr", "magrittr", "stringr", "ggplot2", "ggfortify", "ggthemes", "ggrepel", "stringi", "tidyverse"), repos="https://cran.rstudio.com")'

# Step 3: Install and test Snakemake

RUN conda install wheel &&\
    conda install snakemake

RUN git clone https://github.com/JasonCharamis/QuteS.git &&\
    git clone https://github.com/JasonCharamis/Snakemake_utils.git

ENV PYTHONPATH="$PYTHONPATH:Snakemake_utils"

RUN curl -O -L https://raw.githubusercontent.com/trinityrnaseq/trinityrnaseq/master/Analysis/DifferentialExpression/run_DE_analysis.pl &&\
    curl -O -L https://raw.githubusercontent.com/trinityrnaseq/trinityrnaseq/master/Analysis/DifferentialExpression/analyze_diff_expr.pl

RUN mv *.pl QuteS/workflow/scripts

ENV PATH="~/.local:$PATH"
ENV PATH="QuteS:$PATH"
ENV PATH="QuteS/workflow/scripts:$PATH"

CMD echo "Everything installed!"

