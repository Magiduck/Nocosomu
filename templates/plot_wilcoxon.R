#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(magrittr)
library(ggplot2)



gene_wicoxon_results = fread("$gene_level_wilcoxon")

parquet_file_mutations = "${params.project_folder}/${bed_file.simpleName}/${bed_file.simpleName}.parquet"
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
all_genes_mutation_frequencies <- dbGetQuery(con, glue("SELECT DISTINCT associated_gene, 
                                  count(distinct icgc_mutation_id) 
                                  FROM '{parquet_file_mutations}'
                                  GROUP BY associated_gene")) %>% as.data.table %>% .[order(-`count(DISTINCT icgc_mutation_id)`), ]

gene_wicoxon_results[all_genes_mutation_frequencies, number_of_mutations := `count(DISTINCT icgc_mutation_id)`, on = .NATURAL]

p <- ggplot(gene_wicoxon_results, aes(number_of_mutations, -log10(w_pvalue))) + 
  geom_point(aes(col=is_bf_significant)) +
  scale_color_manual(values = c("#CACACA", "#FFD700")) +
  ylab("-log10 wilcoxon p-value") +
  xlab("Distinct mutations in gene region") +
  guides(color="none") +
  theme_bw()



fwrite(gene_wicoxon_results, "${model_pickle}_${analysis_type}_table.tsv", sep = '\\t')
ggsave(filename = glue("${gene_level_wilcoxon.baseName}_vs_frequency.png") %>% glue, p, width = 7, height = 5)


