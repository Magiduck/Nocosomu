# Nocosomu: importDataset Module

## Description
The `importDataset` module is responsible for ingesting raw somatic variant datasets and harmonizing their genomic coordinates. It systematically removes putative artifacts and germline variants through a sequential quality control pipeline executed via custom Python and R scripts.

## Core Architecture
The module relies on a structured directory of scripts and templates to execute data transformations.
*   **Nextflow Configuration (`nextflow.config`):** Defines parameter defaults, genome build resource paths, and container mappings.
*   **Workflow Logic (`main.nf` & `processes.nf`):** Manages channel routing, sequential filtering processes, and parallelization across defined chromosomes.
*   **Python Utilities (`resources/usr/bin/`):** Contains standalone executable scripts including `remove_below_coverage.py` and `remove_germline_icgc.py` to filter variants based on provided reference databases.
*   **R Templates (`templates/`):** Utilizes DuckDB and data.table within scripts like `elegible_ICGC.R` to parse Parquet datasets and compute regional mutation burden.

## Input Parameters
The module requires specific configuration parameters to define the target data and execution limits.
*   `data_source`: Dictates the parsing logic, accepting either "ICGC" (Parquet format) or "GE" (TSV format).
*   `dataset`: The file path to the raw somatic variant dataset.
*   `bed`: The file path to the target genomic regions.
*   `chromosomes`: A comma-separated list of chromosomes utilized to split and parallelize the data processing.
*   `subgroup_file`: A TSV file defining cancer subtypes and their associated project codes.
*   `minimum_donors`: The minimum required count of distinct donors harboring a mutation for a region to be deemed eligible.
*   `germlineAlleleFrequency`: The maximum acceptable allele frequency threshold used to identify and remove putative germline variants.

## Quality Control Configuration
The module automatically applies multiple sequence-level filters to remove variants that overlap with problematic genomic regions. The `skipQC` parameter accepts a concatenated string of specific keywords to bypass individual filters dynamically.
*   `COV30`: Bypasses the removal of variants located in regions lacking sufficient sequencing depth.
*   `REPEAT`: Bypasses the removal of variants located within known structural repeat elements.
*   `COMPLEXITY`: Bypasses the removal of variants located within low-complexity sequence tracks.
*   `UMAPK36`: Bypasses the removal of variants located in regions with low multi-read mappability scores.
*   `CDS`: Bypasses the removal of variants overlapping established coding sequences.

## Output Formats
The module generates two finalized outputs for downstream mutational burden analysis.
*   **Filtered Variants:** An output collection named `merged_mutations_in_region.tsv` containing the fully processed somatic variants that passed all active quality control checks.
*   **Eligible Regions:** A comma-separated value file named `elegible_regions.csv` containing the final list of genomic regions that meet the minimum donor mutational burden threshold.
