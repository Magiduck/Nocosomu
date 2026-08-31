#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

#get preliminary files
parquet_file = "$parquet_file"
region_name = "${bed_file.simpleName}"

#load ICGC data
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
dbExecute(con,"SET max_memory = '10Gb'")
dbExecute(con,"SET threads = 1")

subgroup_data <- fread("$subgroup_file")
all_subgroups <- subgroup_data[!is.na(my_cancer_type),unique(my_cancer_type)]

gene_transcript_table <- dbGetQuery(con, glue("SELECT DISTINCT gene_affected,transcript_affected FROM '{parquet_file}' WHERE chromosome = '$chromosome'")) %>% as.data.table

all_chromosomes = c("$chromosome")

load_type1 <- function(table_raw){
  setkey(table_raw, chromStart, chromEnd)
  table_raw[,transcript_name := NA]
  return(table_raw)
}

#load region file
bed_regions_raw = load_type1(fread("$bed_file"))
unimapped_genes = unique(bed_regions_raw[, .(`#chrom`, name)]) %>% .[,.N,by=name] %>% .[N==1, name] #ignore genes in more than 1 chromosome

bed_regions = bed_regions_raw[name %in% unimapped_genes]
bed_regions[, min_pos := min(chromStart), by=name]
bed_regions[, max_pos := max(chromEnd), by=name]

update_table_for_region <- function(donor_mutation_file, regions_file, region_label){
  mutation_gene_table = foverlaps(donor_mutation_file, regions_file, by.x=c("chromosome_start", "chromosome_end"), type = "within", nomatch=NULL, which=TRUE)
  donor_mutation_file[mutation_gene_table[,xid], `:=`(associated_gene = regions_file[mutation_gene_table[,yid], name],
                                                      associated_transcript = regions_file[mutation_gene_table[,yid], transcript_name],
                                                      min_pos_region = regions_file[mutation_gene_table[,yid], min_pos],
                                                      max_pos_region = regions_file[mutation_gene_table[,yid], max_pos],
                                                      gene_region = region_label)]
}

for (i in 1:length(all_subgroups)){
  current_subgroup = all_subgroups[i]
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

      update_table_for_region(donor_mutation_file, bed_regions[strand == current_strand & `#chrom` == glue("chr{current_chr}")], region_name)

      writable_file = donor_mutation_file[gene_region != "no_associated_gene_in_strand",.(icgc_mutation_id,
                                                                                          icgc_donor_id,
                                                                                          chromosome,
                                                                                          chromosome_start,
                                                                                          chromosome_end,
                                                                                          strand,
                                                                                          mutation_type,
                                                                                          reference_genome_allele,
                                                                                          mutated_from_allele,
                                                                                          mutated_to_allele,
                                                                                          subgroup,
                                                                                          associated_gene,
                                                                                          gene_region,
                                                                                          min_pos_region,
                                                                                          max_pos_region
                                                                                          )] %>% unique

      setnames(writable_file, "icgc_mutation_id", "mutation_id")
      setnames(writable_file, "icgc_donor_id", "donor_id")
      setnames(writable_file, "associated_gene", "associated_region")
      setnames(writable_file, "gene_region", "region_metadata")

      fwrite(writable_file, "mutations_in_region.tsv", sep = "\t", append = ifelse(current_subgroup != all_subgroups[1], TRUE, FALSE))
    }
  }
}

mutation_file = fread("mutations_in_region.tsv")
mutation_file2 = fread("mutations_in_region.tsv")
mutation_file2[, subgroup := "all_subtypes"]

merged_mutation_file = rbind(mutation_file, mutation_file2) %>% unique %>% as.data.table

# Filter to include only mutations with a recurrence in icgc_mutation_id of < 12, considering only rows where subgroup does not equal "all_subtypes"
#mutation_counts <- merged_mutation_file[subgroup != "all_subtypes", .N, by=mutation_id]
#filtered_mutation_file <- merged_mutation_file[mutation_id %in% mutation_counts[N < 12, mutation_id]]

#fwrite(filtered_mutation_file, "merged_mutations_in_region.tsv", sep = "\t")
fwrite(merged_mutation_file, "merged_mutations_in_region.tsv", sep = "\t")

