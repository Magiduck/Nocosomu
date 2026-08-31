#!/usr/bin/env Rscript 

library(data.table)
library(magrittr)
library(glue)

#Make ism tsvs for prediction

current_chromosome = "${chromosome}"
gene_selection = fread("$gene_selection", header = F) %>% .[,V1]
bed <- fread("$bed_file")

bed[grepl("chr", `#chr`),`#chr` := gsub("chr","",`#chr`)] #Force format
#bed[grepl("\\\\.", name), name := sapply(strsplit(name, "\\\\."), "[", 1)] #Remove version name

all_genes = bed[`#chr` == current_chromosome & name %in% gene_selection, unique(name)]

nucleotides = c("A","C","G","T")
nucleotides2 = c("A","C","G","T")
mutations = CJ(nucleotides,nucleotides2) %>% .[nucleotides < nucleotides2]
mutations[, `#chr` := current_chromosome]

for (current_gene in all_genes){
  print(current_gene)
  
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
  
  fwrite(final_tsv[, .(`#CHROM`, STA, POS, REF, ALT, ID)], glue("{current_gene}_ism.tsv"), sep = '\\t')
}

if (length(all_genes) == 0){
  fwrite(mutations, glue("NO_GENES_{current_chromosome}_ism.tsv"), sep = '\\t')
}

