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

possible_regions <- fread("elegible_regions.csv") %>% .[,unique(associated_region)]
mutations <- fread("mutations.csv")

fwrite(mutations[associated_region %in% possible_regions,], glue("mutations_to_run.csv"), sep = "\\t")


