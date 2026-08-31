#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

#get preliminary files
parquet_file = "$dataset"
parquet_intermediates = "."
region_name = "${bed_file.simpleName}"

#load ICGC data
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")

subgroup_data <- fread("$params.subgroup_file")
all_subgroups <- subgroup_data[, unique(my_cancer_type)]

gene_transcript_table <- dbGetQuery(con, glue("SELECT DISTINCT gene_affected,transcript_affected FROM '{parquet_file}' WHERE chromosome = '$chromosome'")) %>% as.data.table

all_chromosomes = "$chromosome"

load_type1 <- function(table_raw){
  setkey(table_raw, start, stop)
  table_raw[,transcript_name := NA]
#  table_raw[,name := sapply(strsplit(name, "\\\\."), "[", 1)]
  return(table_raw)
}

#load region file
#bed_regions_raw = load_type1(fread("/groups/umcg-fg/tmp02/projects/non-coding-somatic/region_beds/1200_upstream_300_downstream_regions_of_protein_coding_genes.bed.gz"))
bed_regions_raw = load_type1(fread("$bed_file"))
unimapped_genes = unique(bed_regions_raw[, .(`#chr`, name)]) %>% .[,.N,by=name] %>% .[N==1, name] #ignore genes in more than 1 chromosome

bed_regions = bed_regions_raw[name %in% unimapped_genes]
bed_regions[, min_pos := min(start), by=name]
bed_regions[, max_pos := max(stop), by=name]

update_table_for_region <- function(donor_mutation_file, regions_file, region_label){
  mutation_gene_table = foverlaps(donor_mutation_file, regions_file, by.x=c("chromosome_start", "chromosome_end"), type = "within", nomatch=NULL, which=TRUE)
  donor_mutation_file[mutation_gene_table[,xid], `:=`(associated_gene = regions_file[mutation_gene_table[,yid], name],
                                                    associated_transcript = regions_file[mutation_gene_table[,yid], transcript_name],
                                                    min_pos_region = regions_file[mutation_gene_table[,yid], min_pos],
                                                    max_pos_region = regions_file[mutation_gene_table[,yid], max_pos],
                                                    gene_region = region_label)]
}


for (current_subgroup in all_subgroups){
  print(current_subgroup)
  relevant_projects = subgroup_data[my_cancer_type == current_subgroup, unique(project_code) %>% paste0( collapse = "', '")]
  relevant_specimens = subgroup_data[my_cancer_type == current_subgroup & tissue_type == 'primary_cancer', unique(icgc_specimen_id)]
  for (current_chr in all_chromosomes){
    print(current_chr)
    print("Reading")
    donor_mutation_file = dbGetQuery(con,glue("SELECT distinct * exclude(consequence_type, aa_mutation, cds_mutation, gene_affected, transcript_affected) 
                                              from '{parquet_file}'
                                              where chromosome == '{current_chr}'
                                              and project_code IN ('{relevant_projects}')
                                              ")) %>% as.data.table %>% .[icgc_specimen_id %in% relevant_specimens,]
    donor_mutation_file[, chromosome_start := chromosome_start - 1] # overlaps have to be range mode
    donor_mutation_file[, subgroup := current_subgroup] # annotate subgroup
    
    for (current_strand in c("+", "-")){
      print(current_strand)
      donor_mutation_file[,associated_gene := ""]
      donor_mutation_file[,associated_transcript := ""]
      donor_mutation_file[,strand := current_strand]
      donor_mutation_file[,gene_region := "no_associated_gene_in_strand"]
      
      update_table_for_region(donor_mutation_file, bed_regions[strand == current_strand & `#chr` == glue("chr{current_chr}")], region_name)
      #HOTFIX limit to SBS
      write_parquet(donor_mutation_file[gene_region != "no_associated_gene_in_strand" & mutation_type == 'single base substitution',], glue("{parquet_intermediates}/{current_subgroup}_{current_chr}_{current_strand}.parquet"))
    }
  }
}


