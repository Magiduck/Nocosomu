#!/usr/bin/env Rscript 

library(glue)
library(magrittr)
library(data.table)

fold_result = list.files(path = ".", pattern = "*_predictions*")

all_folds_data = lapply(fold_result, fread) %>% do.call("rbind", .)
fwrite(all_folds_data, glue("cadd_merged.tsv"), sep = '\\t')