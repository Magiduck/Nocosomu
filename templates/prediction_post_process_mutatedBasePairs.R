#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

post_process_prediction <- function(file_to_read){
  
  ism_predictions_raw <- fread(file_to_read)
  #ism_predictions_raw[, ID := variant_id]
  ism_predictions_raw[, diff_pred_avg := logSED]
  ism_predictions_raw[, tri_context_ref := "NNN"]
  ism_predictions_raw[, tri_context_alt := "NNN"]
  ###HOTFIX complete reverse changes
  ism_predictions_half1 <- separate(ism_predictions_raw[,.(ID, diff_pred_avg, tri_context_ref, tri_context_alt)], "ID", c("chr","start","end","from","to"), sep = "_")
  ism_predictions_half2 <- separate(ism_predictions_raw[,.(ID, diff_pred_avg, tri_context_ref, tri_context_alt)], "ID", c("chr","start","end","from","to"), sep = "_")
  setnames(ism_predictions_half2, c("from", "to"), c("to", "from"))
  setnames(ism_predictions_half2, c("tri_context_ref", "tri_context_alt"), c("tri_context_alt", "tri_context_ref"))
  ism_predictions_half2[, diff_pred_avg := -diff_pred_avg]
  ism_predictions <- rbind(ism_predictions_half1, ism_predictions_half2)
  ####
  ism_predictions[, ID := paste(chr, start, end, from, to, sep = "_")]
  
  return(ism_predictions)
}

current_feature = strsplit("${prediction_file.baseName}", "_mutatedBasePairs")[[1]][1] #baseName

write_feather(post_process_prediction("$prediction_file"), glue("{current_feature}_mutatedBasePairs_predictions.feather"))


#fix old predictions
#fix_prediction <- function(file_to_read){
#  print(file_to_read)
#  write_feather(post_process_prediction(file_to_read), gsub("tsv.gz", "feather", file_to_read))
#  return()
#}


#for (prediction_file in list.files("./", recursive = F, pattern = "*tsv.gz")){
#  fix_prediction(prediction_file)
#}
