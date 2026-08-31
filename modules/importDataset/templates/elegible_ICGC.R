#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

#get preliminary files
minimum_donors = $minimum_donors

#get real mutation data
project_gene_mutations_all <- fread("$merged_mutations") %>% 
  .[mutation_type == "single base substitution"] %>% #only snvs
  .[max_pos_region - min_pos_region > 56 & max_pos_region - min_pos_region < 1403] # mutations in regions of reasonable size

#remove SNVs that fall within an MNV on same donor(not actually SNVs)
project_gene_mbs <- fread("merged_mutations_in_region.tsv") %>% .[mutation_type == "multiple base substitution (>=2bp and <=200bp)",.(mutation_id, donor_id, chromosome_start, chromosome_end)] %>% unique
setkey(project_gene_mbs, donor_id, chromosome_start, chromosome_end)
project_gene_mutations <- foverlaps(project_gene_mutations_all, project_gene_mbs, by.x=c("donor_id","chromosome_start", "chromosome_end")) %>% 
  .[is.na(mutation_id)] %>% 
  .[, .(number_of_donors = .N), by=.(associated_region, subgroup, chromosome)]

#determine regions that have at least threshold of mutations
project_gene_mutations[,cancer_type := subgroup]
specific_subtype_usable_genes = project_gene_mutations[, .(number_of_donors = sum(number_of_donors)), by=.(associated_region, chromosome, cancer_type)] %>% .[number_of_donors >= minimum_donors,]

combined_usable_genes = specific_subtype_usable_genes[,.(associated_region, chromosome, cancer_type)]

fwrite(combined_usable_genes, "elegible_regions.csv")
