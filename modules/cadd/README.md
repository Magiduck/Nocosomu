### `Cadd_predict` Module Execution

The `Cadd_predict` workflow extracts pre-computed functional consequence scores for somatic single nucleotide variants by querying a local database of Combined Annotation Dependent Depletion predictions.

#### Input Channels
*   **Variant Files (`vcflikeFiles`):** The primary input channel accepts mutation files formatted in a standardized VCF-like tabular structure alongside their corresponding variant identifiers.
*   **Database Path (`caddModelPath_ch`):** A secondary input channel supplies the file path to the comprehensive Parquet database containing the pre-calculated whole-genome CADD scores.

#### Execution Logic
1.  **Chunking:** The workflow utilizes the `delimitedFileSplitter` process to divide the input variant files into smaller parallelizable shards of exactly 1500 lines before applying the `transpose` operator to flatten the resulting channel.
2.  **Score Extraction (`cadd_tsv_like`):** The parallelized shards are passed to the `cadd_scores.R` script, which establishes an in-memory DuckDB connection to efficiently query the local `whole_genome_SNVs.parquet` database file. The script determines the minimum and maximum genomic positions from the input shard to optimize the database query and retrieves the corresponding raw functional scores. Whenever an observed variant changes the allele back to the reference genome base, the script automatically assigns a score of 0 to indicate a benign functional effect.
3.  **Merging (`cadd_merge_folds`):** The workflow groups the processed prediction shards by their original variant identifier and the statically assigned "CADD" group identifier using the `groupTuple` operator. These grouped shards are subsequently consolidated into a single unified file by the `merge_cadd_folds.R` template script using standard data table row binding operations.

#### Output Channels
*   **Merged Predictions (`vcflikeFile_merged`):** The sub-workflow emits a flattened channel containing the original identifier, the CADD group identifier, and the fully merged TSV file containing the extracted prediction scores.
*   **Published Files:** The parent workflow assigns a final filename and publishes the output dataset, although the current execution logic routes this output to a `borzoi_out` directory rather than a CADD-specific destination.
