#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(gmp)
library(tidyr)
library(ggplot2)

shift_bp = as.integer("${params.shift_bp}")

ism_dir = "${params.project_folder}/${bed_file.simpleName}/${dataset.simpleName}/${model_pickle.simpleName}/artifacts/${params.mode}/"
current_chromosome = "$chromosome"
cancer_type = "$cancer_type"
current_mode = "$params.mode"

#get preliminary files
parquet_file_mutations = "${params.project_folder}/${bed_file.simpleName}/${bed_file.simpleName}.parquet"

if (cancer_type != "all_subtypes"){
  extra_selection = glue("and subgroup LIKE '{cancer_type}%'")
} else{
  extra_selection = ""
}

gene_selection <- fread("$gene_selection")
usable_genes <- gene_selection[cancer_type == get("cancer_type", envir = .GlobalEnv) & chromosome == current_chromosome, unique(associated_gene)]

#get real mutation data
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
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
                                                         max_pos_region
                                                 FROM '{parquet_file_mutations}'
                                                 WHERE chromosome = '{current_chromosome}'{extra_selection}")) %>% as.data.table %>% .[associated_gene %in% usable_genes]
setkey(chromosome_mutations, associated_gene)

all_genes = chromosome_mutations[, unique(associated_gene)]

all_wilcox_against_null <- lapply(all_genes, function(current_gene){ 
  
  #real mutations
  current_mutations <- chromosome_mutations[associated_gene == current_gene]
  current_mutations[, dummy_mutation_id := paste(chromosome, chromosome_start, chromosome_end, mutated_from_allele, mutated_to_allele, sep = "_")]
  
  #create nearby mutations
  nearby_mutations_left = copy(current_mutations)
  nearby_mutations_left[, chromosome_start := chromosome_start - shift_bp]
  nearby_mutations_left[, chromosome_end := chromosome_end - shift_bp]
  nearby_mutations_right = copy(current_mutations)
  nearby_mutations_right[, chromosome_start := chromosome_start + shift_bp]
  nearby_mutations_right[, chromosome_end := chromosome_end + shift_bp]
  
  #remove mutations that overflow window &
  #update mutation ids
  nearby_mutations = rbind(nearby_mutations_left[chromosome_start > min_pos_region, ], nearby_mutations_right[chromosome_end < max_pos_region, ])
  nearby_mutations[, dummy_mutation_id := paste(chromosome, chromosome_start, chromosome_end, mutated_from_allele, mutated_to_allele, sep = "_")]

  #check that unobserved allele changes are not observed in other patients
  any_patient_unobserved_nearby_mutations = nearby_mutations[!dummy_mutation_id %in% current_mutations[, dummy_mutation_id]]
  
  #populate mutations with predictions
  #ism_predictions_raw <- fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/SuRE4n_tss_selection_m300_p100_intersection_K562_MSE_default/artifacts/ISM/ENSG00000077935.17_7_ism_predictions.tsv.gz")
  ism_predictions <- read_feather(glue("{ism_dir}/{current_gene}_ism_predictions.feather"))
  ism_predictions[ ,`:=`(chr = NULL, start = NULL, end = NULL, from = NULL, to = NULL)]
  ###
  #
  any_patient_unobserved_nearby_mutations[ism_predictions, on = "dummy_mutation_id==ID", delta_sure := diff_pred_avg]
  current_mutations[ism_predictions, on = "dummy_mutation_id==ID", delta_sure := diff_pred_avg]
  
  #perform wilcoxon test
  plot_data = rbind(current_mutations[,.(dummy_mutation_id, delta_sure)], any_patient_unobserved_nearby_mutations[,.(dummy_mutation_id, delta_sure)])
  plot_data[, label := ifelse(dummy_mutation_id %in% current_mutations[,dummy_mutation_id], "observed_mutation", "comparison_group")]
  if (current_mode == '_abs'){
    plot_data[, ism_prediction := abs(delta_sure)]
    wilcoxon_result <- wilcox.test(plot_data[label == "observed_mutation",ism_prediction], 
                                   plot_data[label == "comparison_group",ism_prediction], 
                                   alternative = "greater")
    
  } else {
    plot_data[, ism_prediction := delta_sure]
    wilcoxon_result <- wilcox.test(plot_data[label == "observed_mutation",ism_prediction], 
                                   plot_data[label == "comparison_group",ism_prediction], 
                                   alternative = "two.sided")
  }
  
  p <- ggplot(plot_data, aes(ism_prediction)) + 
    geom_histogram(aes(fill=label), position = 'dodge') +
    scale_fill_manual(values=c("comparison_group" = "#CACACA", "observed_mutation" = "#FA8072")) +
    ggtitle(glue("{current_gene} pvalue={wilcoxon_result[['p.value']]}"))
  
  ggsave(glue("{current_gene}.pdf"), p)
  plot_data[, gene := current_gene]

  write_path = glue("plotData_${model_pickle.simpleName}_{cancer_type}_{current_chromosome}_nearbyInPromoter_out.tsv")
  fwrite(plot_data, write_path, append = file.exists(write_path), sep = "\\t", col.names = !file.exists(write_path))
  
  return(data.table(associated_gene = current_gene,
                    w_statistic = wilcoxon_result[["statistic"]],
                    w_pvalue = wilcoxon_result[["p.value"]],
                    method = "nearbyInPromoter",
                    cancer_type = cancer_type
  ))
}) %>% do.call("rbind", .)

significance_threshold = gene_selection[,0.05/length(unique(associated_gene))]
all_wilcox_against_null[, is_bf_significant := FALSE]
all_wilcox_against_null[w_pvalue < significance_threshold, is_bf_significant := TRUE]
fwrite(all_wilcox_against_null, glue("${model_pickle.simpleName}_{cancer_type}_nearbyInPromoter.tsv"), sep = "\\t")



