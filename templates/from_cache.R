#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

#to_predict <- fread("/groups/umcg-fg/tmp02/projects/non-coding-somatic/scratch_cgut/nocosomu_iDriver/work/11/b69d765bbc74dc708b512fb42112e9/ENSG00000159208_ism.tsv")
#cached_predictions <- fread("/groups/umcg-fg/tmp02/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/6000_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/SuRE4n_tss_selection_m300_p100_intersection_K562_MSE_default/artifacts/ISM/ENSG00000159208_ism_predictions.tsv.gz")

to_predict <- fread("$tsv_file")
cached_predictions <- fread("${params.from_cache}/${dataset.simpleName}/${model_pickle.simpleName}/artifacts/${params.mode}/${tsv_file.baseName}_predictions.tsv.gz")

merged_predictions = to_predict[cached_predictions, on=.(ID==ID), nomatch=0]

fwrite(merged_predictions[, .(ID, tri_context_ref, tri_context_alt, pred_ref_fwd, pred_alt_fwd, pred_ref_rev, pred_alt_rev)], glue("${tsv_file.baseName}_predictions.tsv.gz"), sep = "\\t")