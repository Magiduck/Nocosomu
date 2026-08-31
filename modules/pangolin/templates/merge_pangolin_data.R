#!/usr/bin/env Rscript 

library(glue)
library(magrittr)
library(data.table)
library(splitstackshape)
library(tidyr)

fold_result = list.files(path = ".", pattern = "*_predictions.tsv")

all_folds_data = lapply(fold_result, function(file_name){
  fread(cmd=glue("cat {file_name} | sed 's/:,/;/g'"), header = T, fill = T )
}) %>% do.call("rbind", .)
setnames(all_folds_data, c("ID","CHR","POS","REF","ALT","Pangolin"))


pangolin_predictions <- cSplit(all_folds_data, "Pangolin", sep=";", direction = "long") %>% 
  separate(., "Pangolin", c("Gene", "Loss", "Gain", "Warning"), sep="\\\\|") %>%
  separate(., "Loss", c("Lossposs", "loss_score"), sep=":") %>%
  separate(., "Gain", c("Gainposs", "gain_score"), sep=":")

pangolin_predictions[, gain_score := abs(as.numeric(gain_score))]
pangolin_predictions[,loss_score := abs(as.numeric(loss_score))]
pangolin_predictions[, max_score := ifelse( gain_score > loss_score, gain_score, loss_score)]
pangolin_predictions[, Maxpos := ifelse( gain_score > loss_score, Gainposs, Lossposs)]

model_loss = pangolin_predictions[,.(ID, loss_score, Lossposs, Gene)]
setnames(model_loss, c("ID","pangolin_prediction", "pangolin_pos", "pangolin_gene"))
model_gain = pangolin_predictions[,.(ID, gain_score, Gainposs, Gene)]
setnames(model_gain, c("ID","pangolin_prediction", "pangolin_pos", "pangolin_gene"))
model_max = pangolin_predictions[,.(ID, max_score, Maxpos, Gene)]
setnames(model_max, c("ID","pangolin_prediction", "pangolin_pos", "pangolin_gene"))

fwrite(model_loss, glue("lossSplice.tsv"), sep="\\t")
fwrite(model_gain, glue("gainSplice.tsv"), sep="\\t")
fwrite(model_max, glue("maxSplice.tsv"), sep="\\t")
