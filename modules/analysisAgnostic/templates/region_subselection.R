#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)
library(ggplot2)

##################

#if they are reasonable size
possible_regions <- fread("elegible_regions.csv")

subselection = fread("$region_subselection", header = FALSE)

if(nrow(subselection) > 0){
  selection = fread("$region_subselection", header = FALSE) %>% .[,V1]
} else {
  selection = possible_regions[, unique(associated_region)] # if no subselection everithing elegeible
}

fwrite(possible_regions[associated_region %in% selection,], glue("elegible_regions_to_run.csv"), sep = "\\t") 




