### `SureResNet_predict` Module Execution

The `SureResNet_predict` workflow evaluates sequence variants using the SuRE-CNN (PARM) sequence-based deep learning model.

#### Input Channels
*   **Variant Files (`vcflikeFiles_ch`):** An input channel providing VCF-like formatted mutation files paired with a variant identifier.
*   **Model Weights (`sureModelPath_ch`):** An input channel providing the `.pth` model weights file.
*   **Configuration Parameters:** The workflow relies on `$params.reference_genome` for sequence context mapping and `$params.sure_dir` to locate the execution script.

#### Execution Logic
1.  **Chunking:** The input variant files are split into smaller shards of 15,000 lines using the `delimitedFileSplitter` process and the `transpose` operator.
2.  **Model Prediction (`sureResNet_tsv_like`):** The parallelized shards are scored by the `predict_sure_from_tsv.py` script running inside a GPU-allocated Singularity container named `sequence_based_models.sif`. The script execution is hardcoded with `--batchsize 1000` and `--L_max 600`.
3.  **Strand Calculation and Merging (`parm_merge`):** The workflow groups the resulting prediction shards by their original identifier and model name using `groupTuple(by:[0,1])`. The `merge_sure_data.R` script binds the chunks and computes the differential sequence effects. It subtracts the reference predictions from the alternate predictions to generate distinct outputs for `forwardStrand.tsv`, `reverseStrand.tsv`, and `averageStrand.tsv`.
4.  **Channel Flattening:** The workflow applies `transpose` to unnest the three output strand files and utilizes a `map` closure to append the strand identifier directly to the group ID string (e.g., modelName_averageStrand).

#### Output Channels
*   **Merged Predictions (`vcflikeFile_merged`):** The workflow emits a final flattened channel containing a tuple with the original identifier, the dynamically generated group ID (incorporating both the model and strand names), and the respective prediction file.
