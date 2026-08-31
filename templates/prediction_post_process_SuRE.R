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
ism_predictions_raw[, diff_pred_avg := (pred_alt_fwd + pred_alt_rev - pred_ref_fwd - pred_ref_rev)/2.0]

fwrite(ism_predictions_raw, glue("{current_feature}_processed_predictions.tsv"))

