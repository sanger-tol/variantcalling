# sanger-tol/variantcalling: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.1.8](https://github.com/sanger-tol/variantcalling/releases/tag/1.1.8)] - Shang Tang (patch 8) - [2025-07-03]

### Enhancements & fixes

- Fix bug in subworkflow deepvariant_caller

## [[1.1.7](https://github.com/sanger-tol/variantcalling/releases/tag/1.1.7)] - Shang Tang (patch 7) - [2025-06-26]

### Enhancements & fixes

- Update Deepvariant to version 1.9.0 (conda-free)
- Add module [DEEPVARIANT_VCFSTATSREPORT](https://github.com/nf-core/modules/tree/master/modules/nf-core/deepvariant/vcfstatsreport) to generate visual report from Deepvariant
- Compress and index concatenated VCFs

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| Deepvariant | 1.6.1       | 1.9.0       |
| HTSlib      |             | 1.21        |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.

## [[1.1.6](https://github.com/sanger-tol/variantcalling/releases/tag/1.1.6)] - Shang Tang (patch 6) - [2025-02-10]

### Enhancements & fixes

- Allow DeepVariant and BCFtools to handle CSI files
- Deal with genomes >4GB in Minimap

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| nextflow   | 22.10.1     | 23.10.1     |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.

## [[1.1.5](https://github.com/sanger-tol/variantcalling/releases/tag/1.1.5)] - Shang Tang (patch 5) - [2025-01-14]

### Enhancements & fixes

- Fix bug in alignment subworkflow

## [[1.1.4](https://github.com/sanger-tol/variantcalling/releases/tag/1.1.4)] - Shang Tang (patch 4) - [2024-12-05]

### Enhancements & fixes

- Module updates and remove Anaconda references

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| bcftools    | 1.17        | 1.20        |
| blastn      | 2.14.1      | 2.15.0      |
| deepvariant | 1.5.0       | 1.6.1       |
| Python      | 3.8.3       | 3.9.1       |
| samtools    | 1.17        | 1.21        |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.

## [[1.1.3](https://github.com/sanger-tol/variantcalling/releases/tag/1.1.3)] - Shang Tang (patch 3) - [2024-05-24]

### Enhancements & fixes

- Fixed the bug in the filtering of multiple PacBio files

## [[1.1.2](https://github.com/sanger-tol/variantcalling/releases/tag/1.1.2)] - Shang Tang (patch 2) - [2024-03-14]

### Enhancements & fixes

- Bug fix when index fai file given for reference fasta file

### Parameters

This release with the following initial parameters:

| Old parameter | New parameter |
| ------------- | ------------- |

> **NB:** Parameter has been **updated** if both old and new parameter information is present. </br> **NB:** Parameter has been **added** if just the new parameter information is present. </br> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.

## [[1.1.1](https://github.com/sanger-tol/variantcalling/releases/tag/1.1.1)] - Shang Tang (patch 1) - [2024-02-02]

### Enhancements & fixes

- Bug fix when reference fasta file name end with .fa or .fa.gz

### Parameters

This release with the following initial parameters:

| Old parameter | New parameter |
| ------------- | ------------- |

> **NB:** Parameter has been **updated** if both old and new parameter information is present. </br> **NB:** Parameter has been **added** if just the new parameter information is present. </br> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.

## [[1.1.0](https://github.com/sanger-tol/variantcalling/releases/tag/1.1.0)] - Shang Tang - [2023-12-20]

### Enhancements & fixes

- Updated the CI procedure to use "sanger-tol" rather than "nf-core" names.
- Renamed Sanger related Github CI test workflows.
- nf-core template was updated from 2.7 to 2.8.
- Removed BAM/CRAM index files from the sample sheets.
- Made fasta index file optional from the inputs.
- Imported PacBio readmapping sub-workflows from [sanger-tol/readmapping pipeline](https://github.com/sanger-tol/readmapping/). Therefore, the pipeline can run on unaligned BAM/CRAM samples now.
- Use VCFtools to calculate per site nucleotide diversity.
- Use VCFtools to calculate heterozygosity.

### Parameters

This release with the following initial parameters:

| Old parameter | New parameter       |
| ------------- | ------------------- |
| --gzi         |                     |
|               | --vector_db         |
|               | --align             |
|               | --include_positions |
|               | --exclude_positions |

> **NB:** Parameter has been **updated** if both old and new parameter information is present. </br> **NB:** Parameter has been **added** if just the new parameter information is present. </br> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference. Only `Docker` or `Singularity` containers are supported, `conda` is not supported.

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| DeepVariant | 1.4.0       | 1.5.0       |
| samtools    | 1.16.1      | 1.17        |
| bcftools    | 1.16.1      | 1.17        |
| python      | 3.11.0      | 3.11.4      |
| vcftools    |             | 0.1.16      |
| blast       |             | 2.14.1+     |
| gunzip:     |             | 1.10        |
| minimap2    |             | 2.24-r1122  |
| awk         |             | 5.1.0       |
| untar       |             | 1.30        |

> **NB:** Dependency has been **updated** if both old and new version information is present. </br> **NB:** Dependency has been **added** if just the new version information is present. </br> **NB:** Dependency has been **removed** if version information isn't present.

## [[1.0.0](https://github.com/sanger-tol/variantcalling/releases/tag/1.0.0)] - Xia Yu - [2023-05-03]

Initial release of sanger-tol/variantcalling, created with the [nf-core](https://nf-co.re/) tools.

### Enhancements & fixes

- Created with [nf-core](https://github/nf-core/tools) template v2.7.2.
- Allows calling variants using DeepVariant for PacBio long read data.
- Significant speed improvements made by splitting the genome before calling variants.
- Outputs both vcf and gvcf formats.

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
