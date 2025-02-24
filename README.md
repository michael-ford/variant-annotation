# Gene Annotation VCF Pipeline 🚀

This is a small pipeline that leverages [Nextflow](https://www.nextflow.io/) to annotates tumor somatic VCFs with gene information and variant effects. 💥

## Overview 🎯

This pipeline performs the following steps:
- **Normalization:** Filters VCFs to include only PASS variants and normalizes them using **bcftools**.
- **Annotation:** Uses **Ensembl VEP** (via the official Docker container for release 112.0) to annotate variants with gene symbols, mutation consequences, and protein changes.
- **Parsing:** Runs a custom Python script to generate summary TSV files:
  - A complete summary of all annotations.
  - A filtered summary of protein-changing variants.

## Getting Started 🚀

### Prerequisites

- [Nextflow](https://www.nextflow.io/) installed.
- Docker or Singularity (configured in Nextflow).
- A local VEP cache directory (set as an environment variable, e.g., `$VEP_CACHEDIR`).

### Usage

Clone the repository and run the pipeline with the following command:

```bash
nextflow run annotate-vcfs.nf --vcfs "vcfs-strelka/*vcf.gz" -with-singularity --vep_cachedir $VEP_CACHEDIR
