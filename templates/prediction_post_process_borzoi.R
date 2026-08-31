#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

current_feature = strsplit("${prediction_file.baseName}", "_ism")[[1]][1] #baseName

ism_predictions_raw <- fread("$prediction_file")
#ism_predictions_raw[, ID := variant_id]
ism_predictions_raw[, diff_pred_avg := hf_logSED]
ism_predictions_raw[, tri_context_ref := "NNN"]
ism_predictions_raw[, tri_context_alt := "NNN"]

fwrite(ism_predictions_raw[gene == current_feature], glue("{current_feature}_processed_predictions.tsv"))

