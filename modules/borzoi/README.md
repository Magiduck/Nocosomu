### `Borzoi_predict` Module Execution

The `Borzoi_predict` workflow evaluates sequence variants using the Borzoi deep learning model, predicting mutational effects across multiple cross-validation folds.

#### Input Channels
*   **Variant Files (`vcflikeFiles`):** The primary input channel accepting VCF-like formatted mutation files.
*   **Configuration Parameters:** The workflow pulls static references from configuration variables, including model paths (`params.borzoi_folds`), reference files (`params.reference_genome`, `params.borzoi_annotation`, `params.borzoi_targets`), and formatting strings (`params.borzoi_columns`).

#### Execution Logic
1.  **Formatting (`borzoi_format`):** The workflow utilizes `format_borzoi.R` to standardize the incoming variant file headers based on the `$params.borzoi_columns` string (defaulting to `#CHROM,STA,POS,REF,ALT,ID`).
2.  **Chunking:** The formatted variant file is split into smaller parallelizable shards using `delimitedFileSplitter`, governed by the `params.chunkSize` parameter (default 1500).
3.  **Model Prediction (`borzoi_tsv_like`):** The shards are processed by `borzoi_sed.py`. This process maps the shards against multiple model folds (e.g., `f0,f1,f2,f3`) concurrently. Execution occurs within a GPU-allocated Singularity container (`sequence_based_models.sif`), producing intermediate HDF5 files (`sed.h5`).
4.  **Data Extraction (`borzoi_extract_sed`):** The workflow extracts human-readable TSV predictions (`*_predictions.tsv`) from the intermediate HDF5 files using `extract_sed_from_h5_gnomAD2.py`.
5.  **Fold Merging (`borzoi_merge_folds`):** The workflow applies Nextflow channel operators (`transpose`, `map`, and `groupTuple`) to collect the extracted chunks and match them by their original identifier and model fold. It executes `merge_borzoi_folds.R` to consolidate these folds, calculating the mean `hf_logSED` and the maximum `hf_nDi` predictions across the dataset.
6.  **Splice Analysis (Conditional):** The root workflow checks the `$params.borzoiMode` parameter. If set to `"splice"`, the merged predictions are routed through the `borzoi_splice_effect` process. This step uses an R script and a conversion table (`params.borzoi_conversionTable`) to compute regional splice disruption metrics (`delta_nDi`, `cdelta_nDi`, `max_ratio_change`) before file publication.

#### Output Channels
*   **Merged Predictions:** The `Borzoi_predict` sub-workflow emits the `vcflikeFile_merged` channel, which contains the final aggregated Borzoi predictions.
*   **Published Files:** The parent workflow publishes the finalized datasets to either a `borzoi_out` or `borzoi_out_splice` directory, depending on the operational mode.
