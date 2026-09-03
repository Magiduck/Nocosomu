### `Models_predict` Module Execution

The `Models_predict` workflow serves as a wrapper to execute sequence-based model predictions. In its current state, it exclusively routes data to the Borzoi model.

#### Input Channels
*   **Variant Files (`vcflikeFiles`):** The sole input channel accepting mutation files.

#### Execution Logic
1.  **Model Routing:** The workflow reassigns the input channel to a local variable and passes it directly to the imported `Borzoi_predict` sub-workflow.
2.  **Identifier Annotation:** Upon completion of the Borzoi predictions, the workflow applies a `map` operator to append the suffix `_borzoi` to the primary identifier of the output tuple.
3.  **Disabled Models:** The logic for `SureResNet_predict` (targeting K562 and HEPG2) exists in the codebase but remains commented out and inactive.

#### Output Channels
*   **Aggregated Results (`all_results`):** The workflow emits a single channel containing the modified Borzoi prediction tuples.
