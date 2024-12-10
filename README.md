# PGS Catalog Benchmarking Analysis

## Overview

The `pgsc_benchmark` tool aims to facilitate the widespread performance benchmarking
of polygenic scores (PGS) curated in the [PGS Catalog](<https://www.pgscatalog.org>).
It has been designed to run directly on the output files produced by the [PGS Catalog
Calculator](<https://pgsc-calc.readthedocs.io/en/latest/>) (`pgsc_calc`), providing a
reproducible method for evaluating score performance. Users are able to run the tool
using their own study cohorts and phenotype data, generating a consistent and directly
comparable set of performance metrics across all desired PGS. The aim is to provide
the human genomics community with a more complete understanding of each score’s predictive
ability, allowing the best performing PGS to be identified across a range of diverse
populations and ancestries.

## Getting Started

### Download

The `pgsc_benchmark` tool is available to download as either a Docker or Singularity image. Please
see the [user manual](<https://github.com/PGScatalog/pgsc_benchmark/blob/main/User-manual.md>)
for download instructions.

- Docker image
- Singularity image

An R implementation is also provided in case any issues are encountered running Docker or
Singularity on your system. Please see the separate R user manual
[here](<https://github.com/PGScatalog/pgsc_benchmark/blob/main/User-manual-R.md>).

- [R script](<https://github.com/PGScatalog/pgsc_benchmark/blob/main/pgsc_benchmark.R>)

### Usage

In brief, four main steps are involved in the benchmarking process:

1. Calculate PGS for your cohort using the PGS Catalog Calculator
2. Prepare a demographic/phenotype file for your cohort
3. Run the benchmarking tool on the outputs of (1) and (2)
4. Post-processing and data interpretation

A full description of each step is provided in the user manual, along with detailed
explanations of the output files produced.

## Credits

Written by Joel T. Gibson\
Cardiovascular Epidemiology Unit\
Department of Public Health and Primary Care\
Victor Phillip Dahdaleh Heart and Lung Research Institute\
University of Cambridge
