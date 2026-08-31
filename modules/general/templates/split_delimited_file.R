#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)
library(ggplot2)

##################

#get real mutation data
big_file <- fread("$big_file")

nlines = abs($nlines)
nrows_big_file = nrow(big_file)
n_pieces = ceiling(nrows_big_file/nlines)



for (piece in 1:n_pieces){
  start = nlines * (piece - 1) + 1
  end = ifelse(start + nlines <= nrows_big_file, start + nlines - 1, nrows_big_file)
  fwrite(big_file[start:end], glue("piece{piece}"), sep = "\\t")
}




