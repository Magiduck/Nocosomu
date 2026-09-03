### `backgroundModels` Module Execution

The `backgroundModels` module generates matched null distributions (background models) for mutational burden testing. It processes observed somatic mutations and constructs hypothetical mutation sets based on specific genomic and contextual constraints.

#### Input Channels
*   **Background Processes:** The `alternativeAlleles`, `samePromoter`, and `otherTumors` processes accept a standardized tuple: `val(ID)`, `path(mutations_file)`, and `path(ism_file)`.
*   **Region Metadata:** The `region_metadata` process accepts `val(ID)`, `path(mutations_file)`, and `path(elegible_regions)`.

#### Execution Logic
The module executes memory-intensive R scripts via a Singularity container (`curzua_workR.sif`) to compute specific null distributions.
1.  **Alternative Alleles (`alternativeAlleles_background`):** Generates all possible unobserved base pair changes (A, C, G, T) at the exact genomic coordinates of the observed mutations. It filters out any synthetically generated mutations that were actually observed in the patient cohort.
2.  **Same Promoter (`samePromoter_background`):** Scans the *in silico* mutagenesis (ISM) file to identify unobserved mutations within the same promoter region that share the exact trinucleotide sequence context of the observed mutations. It processes these sequentially to manage memory consumption.
3.  **Other Tumors (`otherTumors_background`):** Constructs a background distribution using mutations observed in different cancer subtypes (excluding the current donor's subtype). 
4.  **Baseline Extraction (`no_background`):** Isolates single base substitutions from the observed mutations and maps their corresponding trinucleotide contexts without generating synthetic variants.
5.  **Metadata Generation (`region_metadata`):** Computes a comprehensive list of all possible single nucleotide variants within the defined eligible regions to facilitate downstream sequence predictions.

#### Output Channels
*   **Background Data:** The background generation processes emit two channels per model: a background mutations TSV (the null distribution variants) and a donor mapping TSV (linking the background variants to specific patient IDs).
*   **Baseline Data:** The `no_background` process emits `sbs_no_background.tsv`.
*   **Regional SNVs:** The `region_metadata` process emits `snvs_in_region.tsv`.
