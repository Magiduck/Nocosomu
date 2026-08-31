#!/usr/bin/env Rscript 

library(data.table)
library(magrittr)
library(glue)
options(scipen=999)
#stub_loc = "/groups/umcg-fg/tmp02/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/2024-05-03_refTSS_v4_1_human_coordinate_srtdb_overlap_for_cutoff_22_900_upstream_300_downstream_1bpTSS_sorted_nonoverlapping_hg19/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/21/merged_mutations_in_region.tsv"
#all_regions <- fread(stub_loc) 

#Make ism tsvs for prediction

elegible_regions <- fread("$elegible_regions") %>% .[,unique(associated_region)]

all_regions <- fread("$mutations_file") %>% .[associated_region %in% elegible_regions,.(associated_region, min_pos_region, max_pos_region, chromosome)] %>% unique

nucleotides = c("A","C","G","T")
nucleotides2 = c("A","C","G","T")
mutations = CJ(nucleotides,nucleotides2) %>% .[nucleotides != nucleotides2]
mutations[, `#chr` := "dummy"]

all_final_tsv <- lapply(all_regions[, unique(associated_region)], function(current_region){
  print(current_region)
  
  current_bed = all_regions[associated_region == current_region]
  current_bed[, min_pos_region := min_pos_region + 1] #necessary?
  
  position_data = apply(current_bed, 1, function(x){ x[2]:x[3]}) %>% unlist %>% as.integer %>% as.character
  locations <- data.table(`#chr` = "dummy",
                          POS = position_data)
  
  final_tsv <- mutations[locations, on="#chr==#chr", allow.cartesian=TRUE]
  final_tsv[, STA := as.character(as.integer(POS) - 1)]
  setnames(final_tsv, c("REF", "ALT", "#CHROM", "POS", "STA"))

  final_tsv[, associated_region := current_region]
  
  return(final_tsv)
  
}) %>% do.call("rbind", .)

all_final_tsv[all_regions, `#CHROM` := chromosome, on=.(associated_region == associated_region)]

all_final_tsv[, ID := paste(`#CHROM`,as.integer(STA),as.integer(POS),REF,ALT, sep = "_")]

output_tsv <- all_final_tsv[, .(`#CHROM`, STA, POS, REF, ALT, ID)]
setnames(output_tsv, c("chromosome","chromosome_start","chromosome_end","mutated_from_allele","mutated_to_allele","id"))
fwrite(output_tsv, glue("snvs_in_region.tsv"), sep = '\\t')
