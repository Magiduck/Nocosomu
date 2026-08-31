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

#get real mutation data
chromosome_mutations <- fread("$mutations_file") %>% .[mutation_type == 'single_base_substitution']
setkey(chromosome_mutations, associated_region)

fwrite(chromosome_mutations, glue("single_base_substitutions.tsv"), sep = "\\t")


