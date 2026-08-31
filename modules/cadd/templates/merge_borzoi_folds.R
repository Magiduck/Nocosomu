#!/usr/bin/env Rscript 

library(glue)
library(magrittr)
library(data.table)

fold_result = list.files(path = ".", pattern = "*_predictions*")

all_folds_data = lapply(fold_result, fread) %>% do.call("rbind", .)
fwrite(all_folds_data[,.(hf_logSED = mean(hf_logSED), mean(hf_REF), mean(hf_ALT)), by = .(si, ID, gene)], glue("borzoi_merged_${borzoi_folds}.tsv"))