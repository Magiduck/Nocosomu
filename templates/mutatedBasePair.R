#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)
library(ggplot2)

shift_bp = as.integer("${params.shift_bp}")

current_chromosome = "${chromosome}"
gene_selection = fread("$gene_selection", header = F) %>% .[,V1]
bed <- fread("$bed_file")

#get preliminary files
parquet_file_mutations = "${params.project_folder}/${bed_file.simpleName}/${bed_file.simpleName}.parquet"

bed[grepl("chr", `#chr`),`#chr` := gsub("chr","",`#chr`)] #Force format
#bed[grepl("\\\\.", name), name := sapply(strsplit(name, "\\\\."), "[", 1)] #Remove version name

all_genes = bed[`#chr` == current_chromosome & name %in% gene_selection, unique(name)]

nucleotides = c("A","C","G","T")
nucleotides2 = c("A","C","G","T")
mutations = CJ(nucleotides,nucleotides2) %>% .[nucleotides < nucleotides2]
mutations[, `#chr` := current_chromosome]

#get all mutations in window observed
#get real mutation data
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
dbExecute(con, "SET max_memory = '10Gb'")
chromosome_mutations <- dbGetQuery(con, glue("SELECT DISTINCT icgc_mutation_id,
                                                         chromosome, 
                                                         chromosome_start, 
                                                         chromosome_end, 
                                                         mutation_type, 
                                                         gene_region, 
                                                         associated_gene,
                                                         min_pos_region,
                                                         max_pos_region
                                                 FROM '{parquet_file_mutations}'
                                                 WHERE chromosome = '{current_chromosome}'")) %>% as.data.table %>% .[associated_gene %in% gene_selection]
setkey(chromosome_mutations, associated_gene)

for (current_gene in all_genes){
  print(current_gene)
  
  #get all positions observed for that gene
  #get all positions with shift
  mutatedBasePairs = chromosome_mutations[associated_gene == current_gene, unique(chromosome_end)]
  min_pos_region = chromosome_mutations[associated_gene == current_gene, min(min_pos_region)]
  max_pos_region = chromosome_mutations[associated_gene == current_gene, max(max_pos_region)]
  rightShiftedBasePairs = as.integer(mutatedBasePairs) + shift_bp
  leftyShiftedBasePairs = as.integer(mutatedBasePairs) - shift_bp
  allesShiftedBasePairs = unique(c(leftyShiftedBasePairs, rightShiftedBasePairs))
  
  allValidBasePairs = c(mutatedBasePairs, as.character(allesShiftedBasePairs[allesShiftedBasePairs > min_pos_region & allesShiftedBasePairs < max_pos_region]))
  
  #Multiple disjoint regions per gene not supported
  
  current_bed = bed[name == current_gene]
  current_bed[, start := start + 1]
  
  position_data = apply(current_bed, 1, function(x){ x[2]:x[3]}) %>% unlist %>% as.integer %>% as.character
  locations <- data.table(`#chr` = current_chromosome,
                          POS = position_data)
  
  final_tsv <- mutations[locations, on="#chr==#chr", allow.cartesian=TRUE]
  final_tsv[, STA := as.character(as.integer(POS) - 1)]
  setnames(final_tsv, c("REF", "ALT", "#CHROM", "POS", "STA"))
  final_tsv[, ID := paste(`#CHROM`,as.integer(STA),as.integer(POS),REF,ALT, sep = "_")]
  final_tsv[, associated_gene := current_gene]
  
  fwrite(final_tsv[POS %in% allValidBasePairs, .(`#CHROM`, STA, POS, REF, ALT, ID)], glue("{current_gene}_ism.tsv"), sep = '\\t')
}

if (length(all_genes) == 0){
  fwrite(mutations, glue("NO_GENES_{current_chromosome}_ism.tsv"), sep = '\\t')
}

