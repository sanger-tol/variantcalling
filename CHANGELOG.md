# sanger-tol/variantcalling: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.0.0 - [date]

Initial release of sanger-tol/variantcalling, created with the [nf-core](https://nf-co.re/) template.

### Parameters

This release with the following initial parameters:

| Old parameter | New parameter        |
| ------------- | -------------------- |
|               | --input              |
|               | --fasta              |
|               | --fai                |
|               | --gzi                |
|               | --interval           |
|               | --split_fasta_cutoff |

> **NB:** Parameter has been **updated** if both old and new parameter information is present. </br> **NB:** Parameter has been **added** if just the new parameter information is present. </br> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| DeepVariant |             | 1.4.0       |
| samtools    |             | 1.16.1      |
| bcftools    |             | 1.16.1      |
| nextflow    |             | 22.10.6     |
| python      |             | 3.11.0      |
| python      |             | 3.8.3       |
| pigz        |             | 2.3.4       |
| yaml        |             | 6.0         |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.
