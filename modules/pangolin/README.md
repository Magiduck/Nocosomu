### `Pangolin_predict` Module Execution

The `Pangolin_predict` workflow evaluates sequence variants to predict their impact on RNA splicing utilizing the Pangolin deep learning model.

#### Input Channels
*   **Variant Files (`vcflikeFiles_ch`):** The primary input channel, providing either VCF-formatted or TSV-formatted mutation files alongside a unique variant identifier.
*   **Annotation Database (`pangolinAnnotation_ch`):** An input channel supplying the path to the required `.db` annotation file (e.g., GENCODE).
*   **Configuration Parameters:** The workflow pulls static references for the genome build (`params.reference_genome`) and specific execution arguments (`params.pangolin_params`).

#### Execution Logic
1.  **Format Verification (`pangolin_format`):** The workflow evaluates the `$params.fmt` string. If set to `"vcf"`, the input file routes through the `pangolin_format` process, which uses `bcftools` and `awk` to extract standard SNV coordinate columns (CHROM, POS, REF, ALT) and generates a concatenated identifier string. Otherwise, the workflow passes the TSV file directly to the next stage.
2.  **Chunking:** The formatted variant file is split into smaller parallelizable shards of exactly 15,000 lines using `delimitedFileSplitter`, and empty shards are filtered out.
3.  **Model Prediction (`pangolin_tsv_like`):** The shards are processed by the Pangolin execution script. Execution occurs within a Singularity container (`pangolin_git.sif`). A shell command first translates the incoming TSV format to a comma-separated format before invoking the core Pangolin model using the arguments defined in `$params.pangolin_params` (defaulting to `-m False -d 500 --column_ids 'CHR,POS,REF,ALT'`).
4.  **Parsing and Merging (`pangolin_merge`):** The workflow groups the resulting prediction shards by their original identifier and the statically assigned "pangolin" group ID. The `merge_pangolin_data.R` script uses standard data table operations and string splitting (e.g., `cSplit`) to parse the complex delimited Pangolin output string into distinct metrics for splicing loss and splicing gain. The script generates three distinct output files based on these metrics: `lossSplice.tsv`, `gainSplice.tsv`, and `maxSplice.tsv`.
5.  **Channel Flattening:** The workflow applies the `transpose` operator to unnest the three output metric files and utilizes a `map` closure to append the specific splice identifier directly to the group ID string (e.g., pangolin_lossSplice).

#### Output Channels
*   **Merged Predictions (`vcflikeFile_merged`):** The sub-workflow emits a flattened channel containing the original identifier, the dynamically generated group ID (incorporating both the model and the specific splice metric name), and the respective parsed prediction file.
*   **Published Files:** The parent workflow assigns a final filename and publishes the output dataset to a `Pangolin_out` directory.
