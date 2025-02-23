FROM ubuntu:20.04

LABEL maintainer="Michael Ford <michael.kb.ford@gmail.com>"
LABEL version="1.0"
LABEL description="Base container for pipeline: includes bcftools, tabix, and Python (with pysam)"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    wget \
    ca-certificates \
    bcftools \
    tabix \
    python3 \
    python3-pip && \
    python-is-python3 && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip3 install pysam

WORKDIR /data

CMD ["/bin/bash"]
