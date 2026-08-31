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
mutation_file <- fread("$big_file")
colnames(mutation_file) <- strsplit("$borzoi_columns",",")[[1]]

fwrite(mutation_file[,.(`#CHROM`,POS, POS, REF, ALT, ID)], "output_file.tsv")
