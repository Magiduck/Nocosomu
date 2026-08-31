#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

substrLeft <- function(x, n){
  substr(x, 1, n)
}

alpha_missense_parquet <- "/groups/umcg-fg/tmp01/projects/non-coding-somatic/vip/alpha_missense/AlphaMissense_hg19.parquet"

current_feature = strsplit("${tsv_file.baseName}", "_ism_vep")[[1]][1] #baseName
current_feature_noversion = strsplit(current_feature, "\\\\.")[[1]][1] #hotfix

aminoacid_changes <- fread(cmd = "zcat $tsv_file | grep -Ev '^##'") %>% .[Feature == current_feature_noversion]
aminoacid_changes <- separate(aminoacid_changes, Amino_acids, c("from_aa", "to_aa"), "/") %>% as.data.table
aminoacid_changes[is.na(to_aa), to_aa := from_aa]
aminoacid_changes[, protein_variant := paste0(from_aa, Protein_position, to_aa)]

###HOTFIX complete reverse changes
ism_predictions_half1 <- separate(aminoacid_changes[,.(`#Uploaded_variation`, from_aa, to_aa)], `#Uploaded_variation`, c("chr","start","end","from","to"), sep = "_")
ism_predictions_half2 <- separate(aminoacid_changes[,.(`#Uploaded_variation`, from_aa, to_aa)], `#Uploaded_variation`, c("chr","start","end","from","to"), sep = "_")
setnames(ism_predictions_half2, c("from", "to"), c("to", "from"))
setnames(ism_predictions_half2, c("from_aa", "to_aa"), c("to_aa", "from_aa"))
aminoacid_changes_full <- rbind(ism_predictions_half1, ism_predictions_half2)
####
aminoacid_changes_full[, ID := paste(chr, start, end, from, to, sep = "_")]


#load predictions
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
gene_transcript_table <- dbGetQuery(con, glue("SELECT * FROM '{alpha_missense_parquet}' WHERE transcript_id = '{current_feature}'")) %>% as.data.table
gene_transcript_table[, to_aa := substrRight(protein_variant, 1)]
gene_transcript_table[, from_aa := substrLeft(protein_variant, 1)]

###HOTFIX complete reverse changes
ism_predictions_half1 <- gene_transcript_table
ism_predictions_half2 <- copy(gene_transcript_table)
setnames(ism_predictions_half2, c("REF", "ALT"), c("ALT", "REF"))
setnames(ism_predictions_half2, c("from_aa", "to_aa"), c("to_aa", "from_aa"))
ism_predictions_half2[, am_pathogenicity := 0.0]
ism_predictions_half2[, am_class := "to_reference"]
gene_transcript_table_full <- rbind(ism_predictions_half1, ism_predictions_half2)
####
gene_transcript_table_full[, POS_chr := as.character(POS)]

merged = gene_transcript_table_full[aminoacid_changes_full, on=.(POS_chr = end, to_aa = to_aa), nomatch=0]
merged[, Gene := current_feature]
merged[, ID := paste(chr, start, POS, from, to, sep = "_")]
merged[, tri_context_ref := NA]
merged[, tri_context_alt := NA]

final_table <- merged[,.(ID, am_pathogenicity, tri_context_ref, tri_context_alt)]
setnames(final_table, "am_pathogenicity", "diff_pred_avg")


write_feather(final_table, glue("{current_feature}_ism_predictions.feather"))