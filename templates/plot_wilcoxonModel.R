#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(magrittr)
library(ggplot2)
library(splitstackshape)

gene_wicoxon_results_grouped <- fread("$gene_level_wilcoxon")
setkey(gene_wicoxon_results_grouped, associated_gene)

if(nrow(gene_wicoxon_results_grouped) == 0){
  fwrite(gene_wicoxon_results_grouped, "summary_of_significant_genes.tsv", sep = '\\t')
  fwrite(gene_wicoxon_results_grouped, "summary_of_significant_gene_cancertype_hits.tsv", sep = '\\t')
} else{
  cancer_driver_genes <- fread("$cancer_drivers") %>% .[,unique(ensmbl)]
  gene_annotation <- fread("$gene_annotation")
  setkey(gene_annotation, V1)
  
  gene_wicoxon_results_grouped[, is_cancer_driver := ifelse(associated_gene %in% cancer_driver_genes, T, F)]
  gene_wicoxon_results_grouped[gene_annotation, symbol := V7]
  
  fwrite(gene_wicoxon_results_grouped[order(min_pvalue),.(symbol, associated_gene, min_pvalue, max_pvalue, is_cancer_driver, is_global_hit, number_of_cancer_types, cancer_types, analysis_type)], "summary_of_significant_genes.tsv", sep = '\\t')
  
  gene_wicoxon_results_stacked <- cSplit(gene_wicoxon_results_grouped, "cancer_types", sep = ',', direction = 'long', type.convert = F)
  fwrite(gene_wicoxon_results_stacked[, .(analyses_with_hit = .N), by = .(symbol, associated_gene, is_cancer_driver, cancer_types)][order(-analyses_with_hit),], "summary_of_significant_gene_cancertype_hits.tsv", sep = '\\t')
}
