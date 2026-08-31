#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)
library(ggplot2)

#Prepare allele changes
all_nucleotide_changes = CJ(c("A","C","G","T"),c("A","C","G","T")) %>% .[V1 != V2]
setnames(all_nucleotide_changes, c("mutated_from_allele", "mutated_to_fake_allele"))

##################

current_mutations <- fread("$mutations_file") %>% .[mutation_type == 'single base substitution']
current_mutations[, dummy_mutation_id := paste(chromosome, chromosome_start, chromosome_end, mutated_from_allele, mutated_to_allele, sep = "_")]

  #alternative mutations
  all_unobserved_allele_changes = all_nucleotide_changes[current_mutations, on = "mutated_from_allele==mutated_from_allele", allow.cartesian = TRUE, nomatch=0] %>%
    .[mutated_to_fake_allele != mutated_to_allele,]
  all_unobserved_allele_changes[, background_mutation_id := paste(chromosome, chromosome_start, chromosome_end, mutated_from_allele, mutated_to_fake_allele, sep = "_")]
  
  #if a mutation in background was observed then remove it from background only for samples of that cancer subtype
  unobserved_allele_changes = all_unobserved_allele_changes[!current_mutations, on=.(background_mutation_id == dummy_mutation_id, subgroup == subgroup)]

  #experimental: make GC-specific background
#  unobserved_allele_changes = unobserved_allele_changes[
#  (mutated_from_allele == "A" & mutated_to_allele == "C" & mutated_to_fake_allele == "G") |
#  (mutated_from_allele == "A" & mutated_to_allele == "G" & mutated_to_fake_allele == "C") |
#  (mutated_from_allele == "T" & mutated_to_allele == "C" & mutated_to_fake_allele == "G") |
#  (mutated_from_allele == "T" & mutated_to_allele == "G" & mutated_to_fake_allele == "C") |
#  (mutated_from_allele == "C" & mutated_to_allele == "T" & mutated_to_fake_allele == "A") |
#  (mutated_from_allele == "C" & mutated_to_allele == "A" & mutated_to_fake_allele == "T") |
#  (mutated_from_allele == "G" & mutated_to_allele == "T" & mutated_to_fake_allele == "A") |
#  (mutated_from_allele == "G" & mutated_to_allele == "A" & mutated_to_fake_allele == "T")]
  
  #make table that says which background is from which donor
  donor_corrected_background = unobserved_allele_changes %>% 
    .[, .(donor_id, dummy_mutation_id, background_mutation_id, subgroup, associated_region)]
  setnames(donor_corrected_background, c("donor_id", "real_mutation_id","matching_background_mutation_id", "subgroup", "associated_region"))

  #add trinucleotide context info
  ism_file <- fread("$ism_file") %>% unique
  unobserved_allele_changes[ism_file, on=.(chromosome==chromosome, chromosome_start==chromosome_start, 
                                   chromosome_end==chromosome_end,
                                   mutated_from_allele==mutated_from_allele,
                                   mutated_to_fake_allele==mutated_to_allele), `:=` (tri_context_ref=tri_context_ref, 
                                                                                tri_context_alt=tri_context_alt,
                                                                                reference_allele=reference_genome)]
  
  
  unobserved_allele_changes[, mutated_to_allele := mutated_to_fake_allele]
  unobserved_allele_changes[, mutated_to_fake_allele := NULL]
  unobserved_allele_changes[, id := background_mutation_id]
  unobserved_allele_changes[, dummy_mutation_id := NULL]
  
#remove mutations that are in real mutations to avoid double predictions
fwrite(unobserved_allele_changes[!id %in% current_mutations[,unique(dummy_mutation_id)]], glue("alternativeAlleles_background.tsv"), sep = "\\t")
fwrite(unique(donor_corrected_background), glue("alternativeAlleles_donorMapping.tsv"), sep = "\\t")



