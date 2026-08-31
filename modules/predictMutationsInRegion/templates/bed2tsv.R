#!/usr/bin/env Rscript 

library(data.table)
library(magrittr)
library(glue)

#Make ism tsvs for prediction


bed <- fread("$bed_file")

bed[grepl("chr", `#chr`),`#chr` := gsub("chr","",`#chr`)] #Force format
bed[, start := start + 1] # make 1-based
#bed[grepl("\\\\.", name), name := sapply(strsplit(name, "\\\\."), "[", 1)] #Remove version name

nucleotides = c("A","C","G","T")
nucleotides2 = c("A","C","G","T")
mutations = CJ(nucleotides,nucleotides2) %>% .[nucleotides < nucleotides2]
mutations[, `#chr` := "dummy"]

final_tsv <- lapply(1:5, function(bed_row){
  print(bed_row)
  current_bed = bed[bed_row, ]
  position_data = apply(current_bed, 1, function(x){ x[2]:x[3]}) %>% unlist %>% as.integer %>% as.character
  locations <- data.table(`#chr` = "dummy",
                          POS = position_data)
  
  final_tsv <- mutations[locations, on="#chr==#chr", allow.cartesian=TRUE]
  final_tsv[, `#chr` := current_bed[, `#chr`]]
  
  final_tsv[, STA := as.character(as.integer(POS) - 1)]
  setnames(final_tsv, c("REF", "ALT", "#CHROM", "POS", "STA"))
  final_tsv[, ID := paste(`#CHROM`,as.integer(STA),as.integer(POS),REF,ALT, sep = "_")]
  
  final_tsv[, `#CHROM` := current_bed[, `#chr`]]
  
  return(final_tsv)
  
}) %>% do.call("rbind", .)

fwrite(final_tsv[, .(`#CHROM`, STA, POS, REF, ALT, ID)], glue("variants_in_region.tsv"), sep = '\\t')

