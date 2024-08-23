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

The `pgsc_benchmark` tool is currently distributed as both a Docker image _(recommended)_
and an R script. Please see the relevant user manual for the implementation you plan to use:

- Docker image: User-manual-Docker.md (link)
- R script: User-manual-R.md (link)

In brief, four main steps are involved in the benchmarking process:

1. Calculate PGS for your cohort using the PGS Catalog Calculator
2. Prepare a demographic/phenotype file for your cohort
3. Run the benchmarking tool on the outputs of (1) and (2)
4. Post-processing and data interpretation

A full description of each step is provided in the user manuals, along with detailed
explanations of the output files produced.

## Credits

Written by Joel T. Gibson\
Email: [jtg29@cam.ac.uk](mailto:jtg29@cam.ac.uk)\
Cardiovascular Epidemiology Unit\
Department of Public Health and Primary Care\
Victor Phillip Dahdaleh Heart and Lung Research Institute\
University of Cambridge
