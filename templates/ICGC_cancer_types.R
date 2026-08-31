#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

#get preliminary files
parquet_file_mutations = "${params.project_folder}/${bed_file.simpleName}/${bed_file.simpleName}.parquet"
minimum_donors = $params.minimum_donors

#Check which genes are usable
bed <- fread("$bed_file")
#bed[grepl("\\\\.", name), name := sapply(strsplit(name, "\\\\."), "[", 1)] #Remove version name
unimapped_genes = unique(bed[, .(`#chr`, name)]) %>% .[,.N,by=name] %>% .[N==1, name] #ignore genes in more than 1 chromosome

#Determine gene subselection
if("${params.gene_selection}" != "NOT_PROVIDED"){
  gene_subselection = fread("${params.gene_selection}", header = FALSE) %>% .[,V1]
} else{
  gene_subselection = unimapped_genes
}
#get real mutation data
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
project_gene_mutations <- dbGetQuery(con, glue("SELECT DISTINCT associated_gene, subgroup, chromosome, COUNT(DISTINCT icgc_donor_id)
                                                 FROM '{parquet_file_mutations}'
                                                 GROUP BY associated_gene, chromosome, subgroup")) %>% as.data.table %>% .[associated_gene %in% gene_subselection]
setnames(project_gene_mutations, "count(DISTINCT icgc_donor_id)", "number_of_donors")

#all subtypes genes
all_subtypes_usable_genes <- project_gene_mutations[, .(number_of_donors = sum(number_of_donors)), by = .(associated_gene, chromosome)] %>% .[number_of_donors > minimum_donors]
all_subtypes_usable_genes[, cancer_type := "all_subtypes"]

project_gene_mutations[,cancer_type := subgroup]
specific_subtype_usable_genes = project_gene_mutations[, .(number_of_donors = sum(number_of_donors)), by=.(associated_gene, chromosome, cancer_type)] %>% .[number_of_donors > minimum_donors,]

combined_usable_genes = rbind(all_subtypes_usable_genes[,.(associated_gene, chromosome, cancer_type)], specific_subtype_usable_genes[,.(associated_gene, chromosome, cancer_type)])

fwrite(combined_usable_genes, "subtypes_genes.csv")
fwrite(unique(combined_usable_genes[,.(chromosome, cancer_type)]), "subtypes_chromosome.csv")
fwrite(combined_usable_genes[,.(unique(cancer_type))], "subtypes.txt", sep = "\\t", col.names = FALSE)
fwrite(combined_usable_genes[,.(unique(associated_gene))], "genes.txt", sep = "\\t", col.names = FALSE)
