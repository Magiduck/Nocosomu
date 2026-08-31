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

#real trinucleotide changes
current_mutations <- fread("$mutations_file") %>% .[mutation_type == 'single base substitution']
current_mutations[, dummy_mutation_id := paste(chromosome, chromosome_start, chromosome_end, mutated_from_allele, mutated_to_allele, sep = "_")]

#only read relevant part of file
header_ism = fread("all_snvs_with_trinucleotide.tsv", nrows = 2) %>% colnames
fwrite(current_mutations[,.(dummy_mutation_id)], "desired_mutations.txt" )
all_ism <- fread(cmd="grep -w -f desired_mutations.txt all_snvs_with_trinucleotide.tsv")
setnames(all_ism, header_ism)

current_mutations[all_ism, on=.(dummy_mutation_id==id), `:=`(tri_context_ref = tri_context_ref, tri_context_alt = tri_context_alt)]

all_observed_trinucleotide_changes = unique(current_mutations[, .(donor_id, tri_context_ref, tri_context_alt, min_pos_region, max_pos_region, associated_region)])

all_trinucleotides = all_observed_trinucleotide_changes[, unique(tri_context_ref)]
#all_trinucleotides = all_trinucleotides[!is.na(all_trinucleotides)] # TODO: remove, tmp solution

for (trinucleotide in all_trinucleotides){
  
print(trinucleotide)

#one trinucleotide at a time to avoid oom
#read imperfect subset
all_ism <- fread(cmd=glue("grep -w {{trinucleotide}} all_snvs_with_trinucleotide.tsv", , .open = "{{", .close = "}}"))
setnames(all_ism, header_ism)

#unobserved trinucleotide changes that are within the regions in those patients
all_fromReference_matchingTrinucleotide_changes = all_observed_trinucleotide_changes[all_ism[reference_genome == mutated_from_allele], on=.(tri_context_ref == tri_context_ref, 
                                                                                                                                            tri_context_alt == tri_context_alt,
                                                                                                                                            max_pos_region >= chromosome_end,
                                                                                                                                            min_pos_region <= chromosome_end), nomatch=0]
#get overwritten during join and dont have their values anymore
all_fromReference_matchingTrinucleotide_changes[, min_pos_region := NULL]
all_fromReference_matchingTrinucleotide_changes[, max_pos_region := NULL]

all_fromReference_matchingTrinucleotide_changes[,chromosome_end := as.character(as.integer(chromosome_start + 1))]
#remove mutations that are in real mutations to avoid double predictions
fwrite(all_fromReference_matchingTrinucleotide_changes[!id %in% current_mutations[,unique(dummy_mutation_id)]], 
       glue("samePromoter_background.tsv"), sep = "\\t", append = trinucleotide != all_trinucleotides[1])

#starting from the current mutations with full group selection reconstruction full donor table (eg al_subtypes and single subtypes)
donor_corrected_background = unique(current_mutations[,.(donor_id, dummy_mutation_id, tri_context_ref, tri_context_alt, subgroup, associated_region)]) %>% 
  .[all_fromReference_matchingTrinucleotide_changes[,.(id, tri_context_ref, tri_context_alt, associated_region)], on=.(tri_context_ref == tri_context_ref, 
                                                                                                    tri_context_alt == tri_context_alt,
                                                                                                    associated_region == associated_region), nomatch=0, allow.cartesian=TRUE] %>% 
  .[dummy_mutation_id != id, .(donor_id, dummy_mutation_id, id, subgroup, associated_region)] %>% unique()
setnames(donor_corrected_background, c("donor_id", "real_mutation_id","matching_background_mutation_id", "subgroup", "associated_region"))

#if a mutation in background was observed then remove it from background only for samples of that cancer subtype
donor_corrected_background_unseen = donor_corrected_background[!current_mutations, on=.(matching_background_mutation_id == dummy_mutation_id, subgroup == subgroup)]


fwrite(unique(donor_corrected_background_unseen), 
       glue("samePromoter_donorMapping.tsv"), sep = "\\t", append = trinucleotide != all_trinucleotides[1])

}
