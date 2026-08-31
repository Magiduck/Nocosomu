#!/usr/bin/env Rscript 

library(data.table)
library(magrittr)
library(glue)
library(ggplot2)
library(duckdb)

out_dir = "${params.project_folder}"
current_gene = "$gene"
current_model = "${model_pickle.simpleName}"
current_data = "${dataset.simpleName}"
current_regions = "${bed_file.simpleName}"

gene_selection = fread("$gene_selection", header = F) %>% .[,V1]

for (current_gene in gene_selection){
  ism_data <- read_feather(glue("{out_dir}/{current_regions}/{current_data}/{current_model}/artifacts/ISM/{current_gene}_ism_predictions.feather"))
  
  aug_data <- fread(cmd = glue("zcat $gencode_annotation | grep {current_gene}")) %>% .[V3=='start_codon', .(V1,V2,V3,V4,V5,V6,V7)] %>% unique
  gene_data <- fread(cmd = glue("zcat $gencode_annotation | grep {current_gene}")) %>% .[V3=='gene', .(V1,V2,V3,V4,V5,V6,V7)] %>% unique
  
  #hotfix Calculate TSS position
  gene_strand = gene_data[,V7]
  tss_pos = ifelse(gene_strand == '+', gene_data[,V4], gene_data[,V5])
  aug_pos = ifelse(gene_strand == '+', aug_data[,V4], aug_data[,V5])
  ism_data[, translation_relative := as.integer(end) - aug_pos]
  ism_data[, translation_relative := translation_relative * ifelse(gene_strand == '+', 1, -1)]
  
  reference_alleles <- ism_data[, .(start, tri_context_ref)]
  reference_alleles[, reference_alleles := substring(tri_context_ref, 1, 1)]
  setnames(reference_alleles, "start", "end")
  reference_alleles_final <- unique(reference_alleles[, .(end, reference_alleles)])
  
  ism_data_from_ref = ism_data[reference_alleles_final, on=.(end == end, from == reference_alleles), nomatch=0]
  
  ism_data_from_ref[, label_plot := as.character(translation_relative)]
  ism_data_from_ref[order(-abs(diff_pred_avg)), rank := .I]
  ism_data_from_ref[rank > 4, label_plot := ""]
  
  p <- ggplot(ism_data_from_ref, aes(as.integer(end), diff_pred_avg, label = label_plot)) + 
    xlab("Genomic position (bp)") +
    ylab("Differential activity") +
    geom_point(aes(col=to)) +
    geom_text()+
    geom_vline(xintercept = aug_pos, colour="#CACACA") +
    geom_vline(xintercept = tss_pos, colour="#999999") +
    ggtitle(current_gene)+
    theme_bw()+
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  ggsave("ism.pdf", p, width = 20, height = 7)
}
