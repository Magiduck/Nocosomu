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

current_mutations <- fread("$mutations_file")  %>% .[mutation_type == 'single base substitution']
ism_file <- fread("$ism_file")
current_mutations[ism_file, on=.(chromosome==chromosome, chromosome_start==chromosome_start, 
                                 chromosome_end==chromosome_end,
                                 mutated_from_allele==mutated_from_allele,
                                 mutated_to_allele==mutated_to_allele), `:=` (tri_context_ref=tri_context_ref, 
                                                                              tri_context_alt=tri_context_alt,
                                                                              reference_allele=reference_genome)]
current_mutations[, id := paste(chromosome, chromosome_start, chromosome_end, mutated_from_allele, mutated_to_allele, sep = "_")]

fwrite(current_mutations, glue("sbs_no_background.tsv"), sep = "\\t")


