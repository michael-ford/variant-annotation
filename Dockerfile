# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04

LABEL maintainer="Mike <your-email@example.com>"
LABEL version="1.0"
LABEL description="Container with VEP, bcftools, tabix, and Python (pysam)"

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Update packages and install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    wget \
    ca-certificates \
    perl \
    git \
    libdbi-perl \
    libdbd-mysql-perl \
    libxml-simple-perl \
    bcftools \
    tabix \
    python3 \
    python3-pip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install Python package pysam
RUN pip3 install pysam

# Install VEP
RUN mkdir -p /opt/vep && cd /opt/vep && \
    git clone https://github.com/Ensembl/ensembl-vep.git && \
    cd ensembl-vep && \
    # Install VEP in non-interactive mode; cache download is skipped (--NO_CACHE)
    perl INSTALL.pl --AUTO a --NO_CACHE 1

# Set environment variables for VEP
ENV PERL5LIB=/opt/vep/ensembl-vep/src
ENV PATH=/opt/vep/ensembl-vep:$PATH

# Set the working directory (can be adjusted as needed)
WORKDIR /data

# Default command (container will run bash)
CMD ["/bin/bash"]
