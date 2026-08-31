#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

ism_predictions_raw <- fread("raw_predictions.tsv")
#ism_predictions_raw[, ID := variant_id]
ism_predictions_raw[, model_predictions := RawScore]
ism_predictions_raw[, tri_context_ref := "NNN"]
ism_predictions_raw[, tri_context_alt := "NNN"]
ism_predictions_raw[, gene := "unknown"]

fwrite(ism_predictions_raw, glue("processed_predictions.tsv"))

