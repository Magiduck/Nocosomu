#!/usr/bin/env Rscript 

library(glue)
library(magrittr)
library(data.table)

fold_result = list.files(path = ".", pattern = "*_predictions.tsv")

all_folds_data = lapply(fold_result, fread) %>% do.call("rbind", .)
all_folds_data[, model_forward := pred_alt_fwd - pred_ref_fwd]
all_folds_data[, model_reverse := pred_alt_rev - pred_ref_rev]
all_folds_data[, model_average := (pred_alt_fwd + pred_alt_rev - pred_ref_fwd - pred_ref_rev)/2.0]

model_average = all_folds_data[,.(ID, tri_context_ref, tri_context_alt, model_average)]
setnames(model_average, c("ID","tri_context_ref", "tri_context_alt", "model_predictions"))
model_reverse = all_folds_data[,.(ID, tri_context_ref, tri_context_alt, model_reverse)]
setnames(model_reverse, c("ID","tri_context_ref", "tri_context_alt", "model_predictions"))
model_forward = all_folds_data[,.(ID, tri_context_ref, tri_context_alt, model_forward)]
setnames(model_forward, c("ID","tri_context_ref", "tri_context_alt", "model_predictions"))

fwrite(model_average, glue("averageStrand.tsv"))
fwrite(model_reverse, glue("reverseStrand.tsv"))
fwrite(model_forward, glue("forwardStrand.tsv"))