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

model_file = "SuRE_ResNet_Attention_SuRE4n_tss_selection_m300_p100_stranded_EnhA_intersection_intersection_K562_TSS_EnhA_strand_False_trial_20_model_CNN"
mode = "ISM"
#model_file = "borzoi_gtexBlood"
#mode = "mutatedBasePairs"

bed_file = "1200_upstream_300_downstream_regions_of_protein_coding_genes"

project_folder = "/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/"
dataset_file = "all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS"

parquet_file_mutations = glue("/{project_folder}/{bed_file}/{bed_file}.parquet")


significant_genes = fread(glue("/{project_folder}/{bed_file}/{dataset_file}/{model_file}/plots/{mode}/summary_of_significant_genes.tsv"))
significant_genes_long = cSplit(indt = significant_genes, splitCols = "cancer_types", direction = "long", sep=",")
#current_gene = 'ENSG00000111653'
#current_subtype = 'cutaneous_melanoma'



ism_dir = glue("/{project_folder}/{bed_file}/{dataset_file}/{model_file}/artifacts/{mode}/")

near_tss_mutations <- list()

for(i in 1:nrow(significant_genes_long)){
  
  current_gene = significant_genes_long[i, associated_gene]
  current_subtype = significant_genes_long[i, cancer_types]
  current_symbol = significant_genes_long[i, symbol]
  
  if(current_subtype == 'all_subtypes'){
    extra = ""
  } else{
    extra = glue(" AND subgroup = '{current_subtype}'")
  }
  
#get real mutation data
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
current_mutations <- dbGetQuery(con, glue("SELECT DISTINCT icgc_mutation_id,
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
                                                         strand,
                                                 FROM '{parquet_file_mutations}'
                                                 WHERE associated_gene = '{current_gene}'{extra}")) %>% as.data.table
dbDisconnect(con)
current_mutations[, dummy_mutation_id := paste(chromosome, chromosome_start, chromosome_end, mutated_from_allele, mutated_to_allele, sep = "_")]

current_strand = current_mutations[1, strand]
tss_pos = ifelse(current_strand == '+', 
                 current_mutations[1, max_pos_region - 300],
                 current_mutations[1, min_pos_region + 300])

ism_predictions <- read_feather(glue("{ism_dir}/{current_gene}_ism_predictions.feather"))
ism_predictions[ ,`:=`(chr = NULL, start = NULL, end = NULL, from = NULL, to = NULL)]

current_mutations[ism_predictions, on = "dummy_mutation_id==ID", delta_sure := diff_pred_avg]

current_mutations[, chromosome_relative := chromosome_start - tss_pos]

if(current_strand == "-"){
  current_mutations[, chromosome_relative := -chromosome_relative]
}

mutation_counts = current_mutations[, .N, by = icgc_mutation_id]
mutation_counts[current_mutations, on =.(icgc_mutation_id == icgc_mutation_id), chromosome_relative := chromosome_relative]
mutation_counts[current_mutations, on =.(icgc_mutation_id == icgc_mutation_id), delta_sure := delta_sure]

max_change = mutation_counts[,max(abs(delta_sure))]

p <- ggplot(mutation_counts, aes(chromosome_relative, abs(delta_sure))) + geom_point(aes(size=N)) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0) +
  xlim(-6001,300) +
  ylim(0,max_change)+
  xlab("Position relative to TSS") +
  ylab("Pred. expression change (absolute)") +
  scale_size_area(breaks = c(1,5,10,50,100)) +
  labs(size = "Number of donors") +
  ggtitle(glue("{current_subtype} ({current_symbol})")) +
  theme_bw()

near_tss_mutations[[i]] <- current_mutations[chromosome_relative > -1000 & abs(delta_sure) > 0.03]

p
#ggsave(glue("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/6000_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/SuRE_ResNet_Attention_SuRE4n_tss_selection_m300_p100_stranded_EnhA_intersection_intersection_K562_TSS_EnhA_strand_False_trial_20_model_CNN/plots/ISM_abs/gene_plots/{current_gene}_{current_subtype}.png"), 
#       p,
#       width = 5, height = 5)

}

all_near_tss_mutations <- do.call("rbind", near_tss_mutations)

fwrite(all_near_tss_mutations, glue("/{project_folder}/{bed_file}/{dataset_file}/{model_file}/plots/{mode}/donors_with_tss_mutations.tsv"), sep="\t")
                   