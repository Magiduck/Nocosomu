### `analysisAgnostic` Module Execution

The `analysisAgnostic` module coordinates mutational burden testing by comparing sequence-based model predictions of observed somatic variants against matched background null models. 

#### Input Channels
The workflow dynamically constructs its input channels by mapping command-line parameters to expected file paths from upstream modules.
*   **Eligible Regions:** Constructed from `$params.chromosomes`, `$params.outdir`, `$params.bedID`, and `$params.datasetID` to locate `elegible_regions.csv` files.
*   **Somatic Mutations:** Constructed from the same parameters to locate `merged_mutations_in_region.tsv` files.
*   **Reference Constraints:** Utilizes `$params.region_subselection` for region filtering, `$params.reference_genome` for trinucleotide mapping, and `$params.region_aggregation` for output consolidation.

#### Execution Logic
The module enforces a sequential data preparation phase followed by conditional parallel analysis branches.
1.  **Preparation:** The workflow filters the input channels to discard chromosomes lacking mutation data. It applies region and mutation subselection thresholds, extracts local metadata, and annotates every variant with its trinucleotide sequence context using the supplied reference genome.
2.  **Background Generation:** The workflow generates a base distribution of observed mutations (`no_background`). It reads the `$params.analyses` flag to trigger independent background generation branches for `alternativeAlleles`, `samePromoter`, or `otherTumors` background permutations.
3.  **Model Prediction:** The `Models_predict` sub-workflow evaluates all generated sequence variants using the specified deep-learning model architectures. 
4.  **Statistical Testing:** The workflow aggregates the raw predictions and processes them through the `wilcoxon_compare` module. This process computes statistical significance between the observed variant predictions and the background null distributions using jittered permutations. 

#### Output Channels
The workflow publishes three primary data deliverables per evaluated chromosome and model configuration.
*   **`regions_summary.tsv`:** The primary statistical enrichment results utilizing two-sided Wilcoxon tests.
*   **`regions_summary_abs.tsv`:** The absolute effect size summaries utilizing one-sided greater Wilcoxon tests.
*   **`tested_mutations.tsv.gz`:** The raw, variant-level prediction scores and statistical testing labels used during the burden calculation.
