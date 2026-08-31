library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

parquet_file_mutations <- "/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/1200_upstream_300_downstream_regions_of_protein_coding_genes/1200_upstream_300_downstream_regions_of_protein_coding_genes.parquet"

#get real mutation data
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
project_gene_mutations <- dbGetQuery(con, glue("SELECT DISTINCT associated_gene, subgroup, chromosome, COUNT(DISTINCT icgc_donor_id)
                                                 FROM '{parquet_file_mutations}'
                                                 GROUP BY associated_gene, chromosome, subgroup")) %>% as.data.table


setnames(project_gene_mutations, "count(DISTINCT icgc_donor_id)", "number_of_donors")


tumors_of_interest <- c("cutaneous_melanoma", 
                        "prostate_cancer", 
                        "diffuse_large_b_cell_lymphoma", 
                        "breast_cancer_ns", 
                        "follicular_lymphoma", 
                        "colon_cancer_ns")

ggplot(project_gene_mutations[subgroup %in% tumors_of_interest], aes(number_of_donors)) + geom_histogram() + facet_grid(~subgroup)

project_gene_mutations[number_of_donors > 10 & subgroup %in% tumors_of_interest][, .(n_mutated_genes_atleast_10donors = .N), by=subgroup][order(-n_mutated_genes_atleast_10donors)]
project_gene_mutations[number_of_donors > 10][, .(n_mutated_genes_atleast_10donors = .N), by=subgroup][order(-n_mutated_genes_atleast_10donors)]

project_gene_mutations[number_of_donors > 1 & subgroup %in% tumors_of_interest][, .(n_mutated_genes_atleast_1donors = .N), by=subgroup][order(-n_mutated_genes_atleast_1donors)]
project_gene_mutations[number_of_donors > 5 & subgroup %in% tumors_of_interest][, .(n_mutated_genes_atleast_5donors = .N), by=subgroup][order(-n_mutated_genes_atleast_5donors)]
project_gene_mutations[number_of_donors > 10 & subgroup %in% tumors_of_interest][, .(n_mutated_genes_atleast_10donors = .N), by=subgroup][order(-n_mutated_genes_atleast_10donors)]
