#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)
library(ggplot2)

final_tsv <- fread("$mutations_file") %>% .[, .(chromosome, as.character(as.integer(chromosome_start)), as.character(as.integer(chromosome_end)), 
                                               mutated_from_allele, mutated_to_allele,
                                               id)] %>% unique
setnames(final_tsv, c("#CHROM", "STA", "POS", "REF", "ALT", "ID"))

fwrite(final_tsv[, .(`#CHROM`, STA, POS, REF, ALT, ID)], glue("vcf_like.tsv"), sep = '\\t')

