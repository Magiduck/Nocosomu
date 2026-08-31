library(data.table)
library(magrittr)
library(glue)
library(stringr)

what = '19'

#Repeat Masker simple repeats
simple_repeats <- fread(cmd = glue("zcat downloads/hg{what}.fa.out.gz | grep Simple_repeat"))
simple_repeats[, repeat_unit_length := str_length(V10) - 3]
simple_repeats[, estimated_repeats := (V7 - V6)/repeat_unit_length]

candidate_positions = simple_repeats[estimated_repeats >= 5 & repeat_unit_length <=5 ][order(estimated_repeats)]

repeat_masker_sr = simple_repeats[,.(V5,V6,V7, estimated_repeats, repeat_unit_length)]
setnames(repeat_masker_sr, c("chromosome", "start", "end", "estimated_repeats", "repeat_unit_length"))

#Biostarts regex homopolymer
#https://www.biostars.org/p/267241/
#cat ../../models/borzoi/hg38.fa | ./a.out > all_homopolymer38.txt
#cat ../../models/borzoi/hg19.fa | ./a.out > all_homopolymer19.txt
biostars_homopolymer <- fread(glue("all_homopolymer{what}.txt")) %>% .[,1:3]
biostars_homopolymer[, estimated_repeats := 6]
biostars_homopolymer[, repeat_unit_length := 1]
setnames(biostars_homopolymer, c("chromosome", "start", "end", "estimated_repeats", "repeat_unit_length"))

#https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg19.hipstr_reference.bed.gz
####HIPSTR
hipSTR_repeats <- fread(cmd = glue("zcat downloads/hg{what}.hipstr_reference.bed.gz")) %>% .[,.(V1,V2,V3,V5,V4)]
setnames(hipSTR_repeats, c("chromosome", "start", "end", "estimated_repeats", "repeat_unit_length"))


###Harmonized selection
all_candidate_repeats = rbind(hipSTR_repeats, biostars_homopolymer, repeat_masker_sr)
all_candidate_repeats[, to_remove := FALSE ]
all_candidate_repeats[estimated_repeats >= 5 & repeat_unit_length <=5, to_remove := TRUE ]

all_candidate_repeats[, start := start - 5]
all_candidate_repeats[, end := end + 5]
setkey(all_candidate_repeats, chromosome, start, end)

repeat_merge <- foverlaps(all_candidate_repeats, all_candidate_repeats[to_remove==TRUE], which = T)

#all_candidate_repeats[repeat_merge[!is.na(yid),unique(xid)],][to_remove==FALSE]
all_candidate_repeats[repeat_merge[!is.na(yid),unique(xid)], to_remove := TRUE]

repeat_selection = all_candidate_repeats[ to_remove == TRUE]

final_blacklist = repeat_selection[, list(chromosome = chromosome, POS = seq(start, end)), by = 1:nrow(repeat_selection)]

final_blacklist[,fwrite(.SD[,.(POS)],glue("artifacts/hg{what}/{.BY}_harmonizedDietlein.txt.gz"), col.names = F), by=chromosome]


