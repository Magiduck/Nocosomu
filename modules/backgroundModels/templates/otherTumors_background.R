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
  all_othertumor_allele_changes = current_mutations[current_mutations[subgroup != "all_subtypes", .(chromosome, chromosome_end, mutated_from_allele, mutated_to_allele, dummy_mutation_id, subgroup, donor_id, associated_region)], 
                                                    on = .(chromosome = chromosome, chromosome_end = chromosome_end, mutated_from_allele=mutated_from_allele, associated_region = associated_region), allow.cartesian = TRUE, nomatch=0] %>%
  .[subgroup != i.subgroup & donor_id != i.donor_id,]
  all_othertumor_allele_changes[, background_mutation_id := i.dummy_mutation_id]
  
  #make table that says which background is from which donor
  donor_corrected_background = all_othertumor_allele_changes %>% 
    .[, .(donor_id, dummy_mutation_id, background_mutation_id, subgroup, associated_region)]
  setnames(donor_corrected_background, c("donor_id", "real_mutation_id","matching_background_mutation_id", "subgroup", "associated_region"))

  #add trinucleotide context info
  ism_file <- fread("$ism_file") %>% unique
  all_othertumor_allele_changes[ism_file, on=.(chromosome==chromosome, chromosome_start==chromosome_start, 
                                   chromosome_end==chromosome_end,
                                   mutated_from_allele==mutated_from_allele,
                                   i.mutated_to_allele==mutated_to_allele), `:=` (tri_context_ref=tri_context_ref, 
                                                                                tri_context_alt=tri_context_alt,
                                                                                reference_allele=reference_genome)]
  
  #mutated to allele has to go
  
  all_othertumor_allele_changes[, mutated_to_allele := i.mutated_to_allele]
  all_othertumor_allele_changes[, i.mutated_to_allele := NULL]
  all_othertumor_allele_changes[, id := background_mutation_id]
  all_othertumor_allele_changes[, dummy_mutation_id := NULL]
  all_othertumor_allele_changes[, i.dummy_mutation_id := NULL]
  all_othertumor_allele_changes[, i.donor_id := NULL]
  #i.mutation_id
  #i.donor_id
  #
  
#remove mutations that are in real mutations to avoid double predictions
#fwrite(all_othertumor_allele_changes[!id %in% current_mutations[,unique(dummy_mutation_id)]], glue("alternativeAlleles_background.tsv"), sep = "\\t")
#write only one mutation to avoid double prediction
fwrite(all_othertumor_allele_changes[1,], glue("otherTumors_background.tsv"), sep = "\\t")
fwrite(unique(donor_corrected_background), glue("otherTumors_donorMapping.tsv"), sep = "\\t")



