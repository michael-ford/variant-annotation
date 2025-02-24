# Gene Annotation VCF Pipeline

This is a small pipeline that leverages [Nextflow](https://www.nextflow.io/) to annotates tumor somatic VCFs with gene information and variant effects. 💥

## Overview 🎯

This pipeline performs the following steps:
- **Normalization:** Filters VCFs to include only PASS variants and normalizes them using **bcftools**.
- **Annotation:** Uses **Ensembl VEP** (via the official Docker container for release 112.0) to annotate variants with gene symbols, mutation consequences, and protein changes.
- **Parsing:** Runs a custom Python script to generate summary TSV files:
  - A complete summary of all annotations.
  - A filtered summary of protein-changing variants.

## Tools & Technologies 🔧

- **Nextflow:** Workflow orchestration for scalability and reproducibility.
- **Docker / Singularity:** Uses a base Docker image, hosted on Docker Hub, as well as the official Ensembl VEP Docker image.
- **bcftools & tabix:** For variant normalization and VCF indexing.
- **Ensembl VEP:** The variant effect predictor (using the official [docker://ensemblorg/ensembl-vep:release_112.0](https://hub.docker.com/r/ensemblorg/ensembl-vep) image).
- **Python & pysam:** For parsing and summarizing the VEP annotations.


## Getting Started 🚀

### Prerequisites

- [Nextflow](https://www.nextflow.io/) installed.
- Singularity
- A local VEP cache directory (set as an environment variable, e.g., `$VEP_CACHEDIR`).

### Usage

Clone the repository and run the pipeline with the following command:

```bash
nextflow run annotate-vcfs.nf --vcfs <WILDCARD-VCF-PATH> -with-singularity --vep_cachedir $VEP_CACHEDIR
