#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

ism_predictions_raw <- fread("raw_predictions.tsv")
ism_predictions_raw[, gene := "unknown"]

fwrite(ism_predictions_raw, glue("processed_predictions.tsv"))

