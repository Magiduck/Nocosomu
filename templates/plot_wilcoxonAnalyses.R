#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(magrittr)
library(ggplot2)

gene_wicoxon_results = fread("$gene_level_wilcoxon") %>% .[is_bf_significant==T] %>% .[order(w_pvalue)]

global_hits <- gene_wicoxon_results[cancer_type == "all_subtypes", unique(associated_gene)]

gene_wicoxon_results_grouped = gene_wicoxon_results[,.(number_of_cancer_types = .N, min_pvalue = min(w_pvalue), max_pvalue = max(w_pvalue), cancer_types = paste(cancer_type, collapse = ",")),by=associated_gene]
gene_wicoxon_results_grouped[, is_global_hit := associated_gene %in% global_hits]
gene_wicoxon_results_grouped[, analysis_type := "$analysis_type"]
#gene_wicoxon_results_grouped[is_global_hit == T, number_of_cancer_types := number_of_cancer_types - 1]

p <- ggplot(gene_wicoxon_results_grouped, aes(number_of_cancer_types, -log10(min_pvalue))) + 
  geom_point(aes(col=is_global_hit)) +
  scale_color_manual(values = c("#CACACA", "#FFD700")) +
  #facet_wrap(~is_cancer_driver)+
  ylab("-log10 wilcoxon p-value") +
  xlab("Number of cancer types") +
  ggtitle("Significant genes vs number of cancer types") +
  guides(color="none") +
  theme_bw()


fwrite(gene_wicoxon_results_grouped, "${model_pickle}_table.tsv", sep = '\\t')
ggsave(filename = glue("${gene_level_wilcoxon.baseName}_acrossTypes.png") %>% glue, p, width = 7, height = 5)


