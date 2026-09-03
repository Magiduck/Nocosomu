### `Alphagenome_predict` Module Execution

The `Alphagenome_predict` workflow evaluates the functional impact of variants using the AlphaGenome sequence-based model. 

#### Input Channels
*   **Variant Files (`vcflikeFiles`):** The workflow accepts an input channel consisting of VCF-like formatted mutation files paired with a unique variant identifier tuple.

#### Execution Logic
1.  **Formatting (`alphagenome_format`):** The workflow routes the input files through the `format_alphagenome.R` script to harmonize genomic coordinates. This process ensures chromosome names utilize a "chr" prefix, converts 0-based BED coordinates to 1-based VCF `POS` coordinates, and generates a standardized `variant_id` string.
2.  **Chunking:** The formatted TSV files are split into smaller subsets using the `delimitedFileSplitter` process, which determines shard size based on the `params.chunkSize` variable.
3.  **Model Prediction (`ALPHAGENOME_PREDICT`):** The shards are processed in parallel by the `alphagenome_from_vcf.py` script running within an NVIDIA GPU-enabled Singularity container (`alphagenome.sif`). This script enforces a strict offline environment by overriding default Kaggle hub API calls, forcing the system to load model weights and the `GRCh38.p13.genome.fa` reference file directly from a local cluster cache (`/groups/umcg-fg/tmp04/projects/non-coding-somatic/kaggle_cache`). The script scores variants within a 1MB sequence window and filters the output strictly to target colorectal biosamples, such as "colonic mucosa", "HCT116", and "Caco-2". It aggregates these filtered scores into a single `colorectal_average` value for each variant-gene interaction.
4.  **Merging (`MERGE_ALPHAGENOME_CHUNKS`):** The workflow uses `groupTuple` to collect the parallelized CSV prediction chunks by their original identifier. It assigns a static "alphagenome" group ID and utilizes an embedded data.table R script to bind the chunks back into a single consolidated file.

#### Output Channels
*   **Merged Predictions (`predictions`):** The workflow emits a final channel containing a tuple with the original identifier, the "alphagenome" group ID, and the merged CSV file (`${groupID}_merged.csv`).
