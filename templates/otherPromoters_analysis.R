#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)

##################

ism_dir = "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/artifacts/${params.mode}/"
current_chromosome = "$chromosome"
cancer_type = "$cancer_type"
n_permutations = 1000

set.seed(12345)

#get preliminary files
parquet_file_mutations = "${params.project_folder}/${bed_file.simpleName}/${bed_file.simpleName}.parquet"

if (cancer_type != "all_subtypes"){
  extra_selection = glue("and subgroup LIKE '{cancer_type}%'")
} else{
  extra_selection = ""
}

#usable_genes <- fread("$gene_selection") %>% .[cancer_type == cancer_type & chromosome == current_chromosome, unique(associated_gene)]
usable_genes <- "$associated_gene"
sampleable_genes <- fread("$gene_selection") %>% .[, unique(associated_gene)]

#get real mutation data
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
gene_tss_all <- dbGetQuery(con, glue("SELECT DISTINCT chromosome, 
                                                      associated_gene,
                                                      min_pos_region
                                                 FROM '{parquet_file_mutations}'")) %>% as.data.table %>% .[associated_gene %in% sampleable_genes]
setkey(gene_tss_all, associated_gene)

chromosome_mutations <- dbGetQuery(con, glue("SELECT DISTINCT icgc_mutation_id,
                                                         icgc_donor_id, 
                                                         icgc_specimen_id,
                                                         chromosome, 
                                                         chromosome_start, 
                                                         chromosome_end, 
                                                         mutation_type, 
                                                         mutated_from_allele, 
                                                         mutated_to_allele, 
                                                         gene_region, 
                                                         associated_gene,
                                                         min_pos_region,
                                                         max_pos_region,
                                                         strand
                                                 FROM '{parquet_file_mutations}'
                                                 WHERE chromosome = '{current_chromosome}'{extra_selection}")) %>% as.data.table %>% .[associated_gene %in% usable_genes]
setkey(chromosome_mutations, associated_gene)

all_genes = chromosome_mutations[, unique(associated_gene)]

all_gene_reference_pos = unique(chromosome_mutations[, .(associated_gene, strand, min_pos_region, max_pos_region)])
all_gene_reference_pos[, reference_point := ifelse(strand == "+", max_pos_region, min_pos_region)]

all_wilcox_against_null <- lapply(all_genes, function(current_gene){ 

  #real mutations
  current_mutations <- chromosome_mutations[associated_gene == current_gene]
  current_mutations[, dummy_mutation_id := paste(chromosome, chromosome_start, chromosome_end, mutated_from_allele, mutated_to_allele, sep = "_")]
  all_current_mutations_ids = unique(current_mutations[,dummy_mutation_id])

  #get predictions
  read_gene <- function(gene_to_read){
    print(gene_to_read)
    #ism_predictions <- fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/1200_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC/SuRE4n_tss_selection_m300_p100_intersection_K562_MSE_default/artifacts/ISM/ENSG00000164362_ism_predictions.tsv.gz")
    #ism_predictions <- fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/1200_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC/SuRE4n_tss_selection_m300_p100_intersection_K562_MSE_default/artifacts/ISM/ENSG00000164362_ism_predictions.feather")
    ism_predictions <- read_feather(glue("{ism_dir}/{gene_to_read}_ism_predictions.feather"))
    ism_predictions[, associated_gene := gene_to_read]

    return(ism_predictions)
  }
  
  #read predictions
  ism_predictions <- read_gene(current_gene)
  ism_predictions_other_genes <- lapply(sample(sampleable_genes[which(sampleable_genes != current_gene)], n_permutations), read_gene) %>% do.call("rbind", .) %>% unique
  
  #calculate distance relative to tss
  current_mutations[all_gene_reference_pos, on=.(associated_gene == associated_gene), reference_p := reference_point]
  ism_predictions_other_genes[all_gene_reference_pos, on=.(associated_gene == associated_gene), reference_p := reference_point]
  current_mutations[, relative_pos := reference_p - chromosome_start]
  ism_predictions_other_genes[gene_tss_all, on=.(associated_gene == associated_gene), relative_pos := reference_p - as.integer(start)]
  
  #Match trinucleotide context and donor
  current_mutations[ism_predictions, on = .(dummy_mutation_id == ID), tri_context_ref := tri_context_ref] # add trinucleotide context to real mutations
  current_mutations[ism_predictions, on = .(dummy_mutation_id == ID), tri_context_alt := tri_context_alt] # add trinucleotide context to real mutations
  matched_mutations <- current_mutations[,.(icgc_mutation_id, icgc_donor_id, tri_context_ref, tri_context_alt, associated_gene, relative_pos)][ism_predictions_other_genes, on = .(tri_context_ref == tri_context_ref, tri_context_alt == tri_context_alt), nomatch = 0, allow.cartesian=TRUE]
  matched_mutations[, similarity_score := abs(relative_pos - i.relative_pos)]
  
  #Get best position per mutation donor pair in OtherGene
  other_gene_mutations = matched_mutations[, .SD[order(similarity_score)][1,], by = .(icgc_mutation_id, icgc_donor_id, associated_gene, tri_context_ref, tri_context_alt, i.associated_gene)]

  #populate predictions
  current_mutations[ism_predictions, on = "dummy_mutation_id==ID", delta_sure := diff_pred_avg]
  other_gene_mutations[, delta_sure := diff_pred_avg]
  
  #perform wilcoxon test
  wilcoxon_result <- wilcox.test(abs(current_mutations[,delta_sure]), 
                                 abs(other_gene_mutations[,delta_sure]), 
                                 alternative = "greater")
  
  return(data.table(associated_gene = current_gene,
                    w_statistic = wilcoxon_result[["statistic"]],
                    w_pvalue = wilcoxon_result[["p.value"]],
                    method = "otherPromoters",
                    cancer_type = cancer_type
                    ))
}) %>% do.call("rbind", .)

significance_threshold = 0.05/as.numeric("$params.number_of_tests")
all_wilcox_against_null[, is_bf_significant := FALSE]
all_wilcox_against_null[w_pvalue < significance_threshold, is_bf_significant := TRUE]
fwrite(all_wilcox_against_null, glue("${model_pickle.simpleName}_{cancer_type}_otherPromoters.tsv"), sep = "\\t")



