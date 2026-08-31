# Nocosomu: Non-Coding Somatic Mutational Burden Pipeline

## Description

Nocosomu is a Nextflow-based bioinformatics pipeline designed to evaluate the functional impact of non-coding somatic mutations. The pipeline prioritizes genomic regions based on mutational burden by comparing observed single nucleotide variants (SNVs) against matched background null models. It leverages sequence-based deep learning models, including PARM (SuRE-CNN), Borzoi, Enformer, and AlphaGenome, to predict differential gene expression effects resulting from these non-coding variants. For coding variants, the pipeline integrates CADD predictions.

## Prerequisites

- Nextflow: (version 22.10.1 or later recommended).

- Singularity/Apptainer: Required for reproducible execution via containerized environments.

- Compute Infrastructure: The pipeline is optimized for SLURM workload managers with GPU availability (NVIDIA A40/V100 recommended for Borzoi and AlphaGenome).

## Methodology

Nocosomu applies a rigorous quality control framework to limit false somatic mutation calls. It anchors promoter regions around refTSS v4.1 transcription start sites (900-bp upstream, 500-bp downstream). Variants are heavily filtered by removing those that overlap with:

- ENCODE blacklists

- gnomAD low-coverage regions

- gnomAD putative germline variants (FAF > 0.01%)

- Simple repeats (RepeatMasker) and short tandem repeats (HipSTR)

- Coding sequences (unless explicitly running the coding pipeline)

To evaluate mutational burden, the pipeline constructs two distinct background null distributions:

1. Alternative Alleles: Compares observed mutations against unobserved alternative base pair changes at the exact same position (e.g., to capture gain-of-function motif creations).

2. Same Promoter: Evaluates base pair changes elsewhere in the same promoter sharing the identical trinucleotide context (e.g., to capture loss-of-function motif disruptions).

Statistical significance is determined using a Wilcoxon rank-sum test, strictly calibrated with 100,000 permutations per region.

## Pipeline Architecture

The Nocosomu codebase executes via standalone module entry points rather than a single root script, allowing for modular dataset processing and analysis.

### Core Modules:

- importDataset: Ingests variant datasets (e.g., Genomics England, HMF, ICGC) and performs strict quality control filtering using the specified genome build resources.

- analysisAgnostic / vcfPrediction: Coordinates mutational burden testing, null model generation, and hypothesis testing.

- backgroundModels: Generates the "Alternative Alleles" and "Same Promoter" mutational backgrounds.

- Model Modules: borzoi, sureResnet (PARM), alphagenome, cadd orchestrate the sequence-based model inference.

## Configuration Profiles

The pipeline relies on several configuration profiles to manage environments and genome builds. You can chain these together using the -profile flag.

- Genome Builds (hg19, hg38): Crucial profiles that automatically map the paths to the reference genomes, ENCODE blacklists, gnomAD datasets, AnySTR/HipSTR resources, and UMAP mappability tracks required for the QC steps.

- Execution Engines (slurm, local_execution): Defines task resources, queues, and retry logic.

- Containers (sif): Mounts the required Singularity images (e.g., curzua_workR.sif for R scripts, sequence_based_models.sif / borzoi.sif / alphagenome.sif for deep learning predictions).

- Institutional (umcg, nibbler, gearshift): Specific cluster configurations for the UMCG environment containing hardcoded reference paths and QoS settings.

## Quick Start / Execution Instructions

### 1. Dataset Import and QC

This command processes and filters the cancer variant datasets.

````
nextflow run modules/importDataset/main.nf \
-profile default_params,umcg,sif,hg38,slurm \
--dataset /groups/umcg-fg/tmp02/projects/non-coding-somatic/cancer/genomics_england/hmf_regions_frequency_matrices_QCed/gel_all_pseudodonorQCed.tsv.gz \
--bed /groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/reftss/artifacts/2024-05-03_refTSS_v4_1_human_coordinate_srtdb_overlap_for_cutoff_22_900_upstream_300_downstream_1bpTSS_sorted_nonoverlapping.bed \
--outdir /groups/umcg-fg/tmp02/projects/non-coding-somatic/umcg-tvanlieshout/pipeline_output/ \
--minimum_donors 5 \
--data_source GE \
--chromosomes "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"
```

### 2. Agnostic Mutational Burden Analysis

This command initiates the sequence-based predictions and mutational burden computations against the background models.

```
nextflow run modules/analysisAgnostic/main.nf \
-profile slurm,default_params,umcg,sif,hg38 \
--datasetID gel_all_pseudodonorQCed \
--bedID 2024-05-03_refTSS_v4_1_human_coordinate_srtdb_overlap_for_cutoff_22_900_upstream_300_downstream_1bpTSS_sorted_nonoverlapping \
--outdir /groups/umcg-fg/tmp02/projects/non-coding-somatic/umcg-tvanlieshout/pipeline_output/ \
--model PARM-K562,PARM-HEPG2 \
--analyses alternativeAlleles,samePromoter \
--reference_genome /groups/umcg-fg/tmp02/projects/non-coding-somatic/models/sure/surecnn/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--chromosomes "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" \
-resume

```
