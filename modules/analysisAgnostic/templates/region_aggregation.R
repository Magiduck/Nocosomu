#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

#concatenate all files
background_mutations = lapply(list.files(path = ".", pattern = "background_mutations*"), fread) %>% do.call("rbind", .)
real_mutations = lapply(list.files(path = ".", pattern = "real_mutations*"), fread) %>% do.call("rbind", .)
background = lapply(list.files(path = ".", pattern = "^background[[:digit:]]"), fread) %>% do.call("rbind", .)
no_background = lapply(list.files(path = ".", pattern = "no_background*"), fread) %>% do.call("rbind", .)
region_selection = lapply(list.files(path = ".", pattern = "region_selection*"), fread) %>% do.call("rbind", .)

#substitute associated reigon id with agreggated id NAs are removed before writing
region_aggregation = fread("$region_aggregation", header = F)

background[region_aggregation, on=.(associated_region==V1), aggregated_region := V2]
no_background[region_aggregation, on=.(associated_region==V1), aggregated_region := V2]

background[, associated_region := aggregated_region]
no_background[, associated_region := aggregated_region]

background[, aggregated_region := NULL]
no_background[, aggregated_region := NULL]


region_selection_union = region_selection[region_aggregation, on=.(associated_region==V1)]
#even if aggregation covers 1 element it goes thorugh (can be changed later)
region_selection_updated = region_selection_union[,.N,by=.(cancer_type,V2)][!is.na(cancer_type),]
region_selection_updated[, chromosome := "dummy"]

region_selection_final = region_selection_updated[,.(V2, chromosome, cancer_type)]
setnames(region_selection_final, "V2", "associated_region")

fwrite(background_mutations, "background_mutations.tsv", sep = "\\t")
fwrite(real_mutations, "real_mutations.tsv", sep = "\\t")
fwrite(background, "background.tsv", sep = "\\t")
fwrite(no_background, "no_background.tsv", sep = "\\t")

fwrite(region_selection_final, "region_selection.tsv", sep = "\\t")