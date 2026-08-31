#!/usr/bin/env Rscript 
library(glue)
library(magrittr)
library(data.table)

#m <- fread("/groups/umcg-fg/tmp02/projects/non-coding-somatic/umcg-curzua/scratch_cgut/borzoi_predict/work/ef/5dbc9e2f5bec7057b025331e58f27f/borzoi_merged_f0.tsv")
m <- fread("$prediction_file") %>% unique

#conversion_table = fread("/groups/umcg-fg/tmp02/projects/non-coding-somatic/shared/gencode/artifacts/gencode.v45lift37.intronsToExonsAndTranscript.tsv.gz")
conversion_table = fread("$conversion_table")
#remove version and create intron_id
#conversion_table[grepl("\\.", exon_5end), exon_5end := sapply(strsplit(exon_5end, "\\."), "[", 1)]
#conversion_table[grepl("\\.", exon_3end), exon_3end := sapply(strsplit(exon_3end, "\\."), "[", 1)]
conversion_table[grepl("\\\\.", exon_5end), exon_5end := sapply(strsplit(exon_5end, "\\\\."), "[", 1)]
conversion_table[grepl("\\\\.", exon_3end), exon_3end := sapply(strsplit(exon_3end, "\\\\."), "[", 1)]
conversion_table[, intron_id := gsub("ENSE","ENSI", exon_5end)]

#annotate m
#annotate transcript
m[, feature_type := ifelse(grepl("ENSE", gene), "exon", "intron")]
m[conversion_table, on=.(gene==exon_5end), transcript := transcript_id]
m[conversion_table, on=.(gene==exon_3end), transcript := transcript_id]
m[conversion_table, on=.(gene==intron_id), transcript := transcript_id]  #intron dummy id is its 5end exon
#specify intron features
m[conversion_table, on=.(gene=intron_id), feature_5end := exon_5end]
m[conversion_table, on=.(gene=intron_id), feature_3end := exon_3end]
#specify exon features
m[conversion_table, on=.(gene=exon_5end), feature_3end := intron_id]
m[conversion_table, on=.(gene=exon_3end), feature_5end := intron_id]
#neighboring exon
m[conversion_table, on=.(gene=exon_3end), previous_exon := exon_5end]
m[conversion_table, on=.(gene=exon_5end), next_exon := exon_3end]

m[m, on=.(si, ID, next_exon==gene, feature_type), ratio_to_next := abs(hf_logSED)/abs(i.hf_logSED)]
m[m, on=.(si, ID, previous_exon==gene, feature_type), ratio_to_previous := abs(hf_logSED)/abs(i.hf_logSED)]

m[is.na(ratio_to_next), max_ratio_change := ratio_to_previous]
m[is.na(ratio_to_previous), max_ratio_change := ratio_to_next]
m[!is.na(ratio_to_next) & !is.na(ratio_to_previous), max_ratio_change := ifelse(ratio_to_previous > ratio_to_next, ratio_to_previous, ratio_to_next)]

#merged_comparison = m[m, on=.(transcript==transcript, ID==ID), allow.cartesian=T][gene != i.gene] %>% 
#  .[, .(local_nDi = mean(i.hf_nDi)), by=.(si, ID, gene, hf_logSED, hf_nDi, transcript)]

#Compared to other similar featues in same transcript
sum_nDi = m[,.(sum_nDi = sum(hf_nDi), nfeatures = .N), by=.(si, ID, transcript, feature_type)]
m[sum_nDi, on=.(si==si, ID==ID, transcript==transcript, feature_type = feature_type), `:=` (sum_nDi = sum_nDi, nfeatures = nfeatures)]
m[, local_nDi := (sum_nDi - hf_nDi)/nfeatures]
m[, delta_nDi := hf_nDi - local_nDi]
m[, cdelta_nDi := hf_nDi/local_nDi]

#finish printing this 
fwrite(m[,.(si, ID, gene, hf_logSED, hf_nDi, V3, V4, transcript, delta_nDi, cdelta_nDi, max_ratio_change)], "spliceEffect_$prediction_file", sep="\\t")

