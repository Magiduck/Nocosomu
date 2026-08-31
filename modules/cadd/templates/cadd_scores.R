#!/usr/bin/env Rscript 
library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
dbExecute(con,"SET max_memory = '10Gb'")
dbExecute(con,"SET threads = 1")

m <- fread("piece.csv")

current_chr = as.character(m[,unique(`#CHROM`)])
min_pos = m[,min(POS)]
max_pos = m[,max(POS)]

results = dbGetQuery(con, paste0('SELECT * exclude(REF) from "whole_genome_SNVs.parquet" as a NATURAL JOIN', ' read_csv("piece.csv", delim=\'\\t\', header = true, columns = {\'#CHROM\': \'VARCHAR\', \'STA\': \'BIGINT\', \'POS\': \'BIGINT\', \'REF\': \'VARCHAR\', \'ALT\':\'VARCHAR\', \'ID\':\'VARCHAR\'}), WHERE a."#CHROM" = ', glue('\'{current_chr}\' and a.POS > {min_pos} and a.POS < {max_pos}'))) %>% as.data.table

#Patients with a change in allele to the reference genome get 0 score (Towards Benign) as compromise
m[results, RawScore := RawScore, on=.(ID==ID)]
m[is.na(RawScore), RawScore := 0]

fwrite(m[,.(`#CHROM`, POS, REF, ALT, ID, RawScore)], "cadd_prediction.tsv", sep = '\\t')

##Chrom      Pos Ref Alt  RawScore  PHRED      STA                      ID
#1       2 25450976   G   C  0.466821  7.540 25450975 2_25450975_25450976_G_C