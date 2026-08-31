#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(magrittr)
library(ggplot2)
library(splitstackshape)

# make plots showcasing the frequencies of ICGC mutations per tumor type
# How many mutations per type

parquet_dir = "/groups/umcg-fg/tmp01/projects/non-coding-somatic/ICGC_per_tumor_type/data_views/columnar_view/"

parquet_expression <- glue("{parquet_dir}/_all_exp_array_ICGC_release28.parquet")
parquet_copy_number <- glue("{parquet_dir}/_all_copy_number_ICGC_release28.parquet")
parquet_consequence <- glue("{parquet_dir}/all_simple_somatic_mutations_ICGC_release28_unsorted.parquet")
#parquet_regions <- "/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/6000_upstream_300_downstream_regions_of_protein_coding_genes/6000_upstream_300_downstream_regions_of_protein_coding_genes_OLD.parquet"

#significant_tss_mutations = fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/1200_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/borzoi_gtexBlood/plots/mutatedBasePairs/donors_with_tss_mutations.tsv")
significant_tss_mutations = fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/6000_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/SuRE_ResNet_Attention_SuRE4n_tss_selection_m300_p100_stranded_EnhA_intersection_intersection_K562_TSS_EnhA_strand_False_trial_20_model_CNN/plots/ISM_abs/donors_with_tss_mutations.tsv")
significant_tss_mutations = fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver//1200_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/SuRE_ResNet_Attention_SuRE4n_tss_selection_m300_p100_stranded_EnhA_intersection_intersection_K562_TSS_EnhA_strand_False_trial_20_model_CNN/plots/ISM/donors_with_tss_mutations.tsv")
setkey(significant_tss_mutations, chromosome, chromosome_start, chromosome_end)
significant_tss_mutations[, chromosome_str := as.character(chromosome)]
significant_tss_donors = paste0("('", significant_tss_mutations[,paste0(unique(icgc_donor_id), collapse = "','")], "')")

#Are genes expressed (Both in microarray (1) and in seq (2) too litle patients overlap with the patients with tss mutations)
#https://docs.icgc.org/dictionary/viewer/#?q=gene_model&viewMode=details&dataType=exp_seq_p
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
gene_expression_donors <- current_mutations <- dbGetQuery(con, glue("SELECT DISTINCT icgc_donor_id, 
                                                 project_code, 
                                                 gene_model,
                                                 gene_id,
                                                 normalized_expression_value,
                                                 platform,
                                                 FROM '{parquet_expression}'
                                                 WHERE icgc_donor_id IN {significant_tss_donors}")) %>% as.data.table
dbDisconnect(con)
gene_expression_donors[, has_tss_mutation := FALSE]
gene_expression_donors[significant_tss_mutations, on = .(gene_affected == associated_gene, icgc_donor_id == icgc_donor_id), has_tss_mutation := TRUE ]


#Are genes intact? (Support for)
non_coding <- "('intergenic_region','intron_variant','upstream_gene_variant','exon_variant','downstream_gene_variant','synonymous_variant')"
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
gene_consequence_donors <- current_mutations <- dbGetQuery(con, glue("SELECT DISTINCT icgc_mutation_id, 
                                                 icgc_donor_id,
                                                 project_code, 
                                                 consequence_type,
                                                 aa_mutation,
                                                 cds_mutation,
                                                 gene_affected,
                                                 transcript_affected
                                                 FROM '{parquet_consequence}'
                                                 WHERE icgc_donor_id IN {significant_tss_donors}
                                                 AND consequence_type NOT IN {non_coding}")) %>% as.data.table
dbDisconnect(con)

gene_consequence_donors[, has_tss_mutation := FALSE]
gene_consequence_donors[significant_tss_mutations, on = .(gene_affected == associated_gene, icgc_donor_id == icgc_donor_id), has_tss_mutation := TRUE ]

projects_from_donors = gene_consequence_donors[has_tss_mutation == T, unique(project_code)]
gene_consequence_donors_subset = gene_consequence_donors[project_code %in% projects_from_donors]

real_number = nrow(gene_consequence_donors[has_tss_mutation == TRUE,])

perm_results <- lapply(1:1000, function(x){
  perm_measurements = gene_consequence_donors_subset[sample(1:nrow(gene_consequence_donors_subset), real_number, replace = F), .N, by = consequence_type]
  perm_measurements[, n_fold := x]
  return(perm_measurements)
  }) %>% do.call("rbind", .)

real_results = gene_consequence_donors[has_tss_mutation == TRUE, .N, by = consequence_type]
real_results[, n_fold := "real"]

perm_plot_data <- rbind(real_results, perm_results)

ggplot(perm_plot_data, aes(consequence_type, log10(N))) + geom_point(aes(col=n_fold == "real")) + coord_flip()

#Is mutation amplified?
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
gene_copynumber_donors <- current_mutations <- dbGetQuery(con, glue("SELECT DISTINCT icgc_donor_id, 
                                                 project_code, 
                                                 copy_number,
                                                 mutation_type,
                                                 segment_mean,
                                                 segment_median,
                                                 chromosome,
                                                 chromosome_start,
                                                 chromosome_end,
                                                 FROM '{parquet_copy_number}'
                                                 WHERE icgc_donor_id IN {significant_tss_donors}")) %>% as.data.table %>% .[!is.na(copy_number)]
gene_copynumber_donors[, chromosome := as.integer(chromosome)]
dbDisconnect(con)
setkey(gene_copynumber_donors, chromosome_str, chromosome_start, chromosome_end)

matched_copy_number = foverlaps(gene_copynumber_donors, significant_tss_mutations, 
          by.x=c("chromosome", "chromosome_start", "chromosome_end"), 
          by.y=c("chromosome", "chromosome_start", "chromosome_end"),
          nomatch=NULL) %>% .[icgc_donor_id == i.icgc_donor_id]

gene_copynumber_donors[,summary(copy_number)]
matched_copy_number[,summary(copy_number)]

gene_copynumber_donors[, has_tss_mutation := FALSE]
gene_copynumber_donors[matched_copy_number, on = .(icgc_donor_id == i.icgc_donor_id, 
                                                   chromosome == chromosome, 
                                                   chromosome_start == i.chromosome_start, 
                                                   chromosome_end == i.chromosome_end), has_tss_mutation := TRUE ]

perm_results_cn <- lapply(1:1000, function(x){
  perm_measurements = gene_copynumber_donors[sample(1:nrow(gene_copynumber_donors), 79, replace = F), .N, by = mutation_type]
  perm_measurements[, n_fold := x]
  return(perm_measurements)
}) %>% do.call("rbind", .)

real_results_cn = matched_copy_number[, .N, by = i.mutation_type]
real_results_cn[, n_fold := "real"]
setnames(real_results_cn, "i.mutation_type", "mutation_type")

perm_plot_data_cn <- rbind(real_results_cn, perm_results_cn)

ggplot(perm_plot_data_cn, aes(mutation_type, N)) + geom_point(aes(col=n_fold == "real")) + coord_flip()# + facet_wrap(~n_fold == "real")

#sure prediction vs copy number
ggplot(matched_copy_number, aes(i.mutation_type, delta_sure)) + geom_violin() + geom_point(position = "jitter")

