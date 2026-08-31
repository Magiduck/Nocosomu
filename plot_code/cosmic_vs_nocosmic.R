#!/usr/bin/env Rscript 

library(duckdb)
library(glue)
library(magrittr)
library(data.table)
library(arrow)
library(magrittr)
library(ggplot2)
library(splitstackshape)

sure_predictions <- fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/1200_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/SuRE_ResNet_Attention_SuRE4n_tss_selection_m300_p100_stranded_EnhA_intersection_intersection_K562_TSS_EnhA_strand_False_trial_20_model_CNN/plots/ISM/alternativeAlleles/all_subtypes/ISM_alternativeAlleles_all_subtypes_vs_frequency.tsv")
all_cosmic_drivers_ensg <- fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/shared/gene_sets/master_driver_sets.csv") %>% .[dataset == "COSMIC CGC",unique(ensmbl)]


wilcox_result = wilcox.test(sure_predictions[associated_gene %in% all_cancer_drivers_ensg, -log10(w_pvalue)],
            sure_predictions[!associated_gene %in% all_cancer_drivers_ensg, -log10(w_pvalue)],
            alternative = "two.sided")
ggplot(sure_predictions, aes(associated_gene %in% all_cosmic_drivers_ensg, -log10(w_pvalue))) + 
  geom_violin(aes(col=associated_gene %in% all_cosmic_drivers_ensg)) +
  geom_boxplot(aes(col=associated_gene %in% all_cosmic_drivers_ensg)) +
  ggtitle(glue("AltAlleles COSMIC vs not COSMIC : wilcox pvalue:{round(wilcox_result$p.value, 3)}")) +
  guides(colour = FALSE) +
  theme_bw()

sure_predictions[, median(-log10(w_pvalue)), by=associated_gene %in% all_cancer_drivers_ensg]
sure_predictions[, mean(-log10(w_pvalue)), by=associated_gene %in% all_cancer_drivers_ensg]

sure_predictions <- fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/1200_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/SuRE_ResNet_Attention_SuRE4n_tss_selection_m300_p100_stranded_EnhA_intersection_intersection_K562_TSS_EnhA_strand_False_trial_20_model_CNN/plots/ISM/nearbyInPromoter/all_subtypes/ISM_nearbyInPromoter_all_subtypes_vs_frequency.tsv")
fwrite(sure_predictions[w_pvalue < 0.05], "/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/1200_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/SuRE_ResNet_Attention_SuRE4n_tss_selection_m300_p100_stranded_EnhA_intersection_intersection_K562_TSS_EnhA_strand_False_trial_20_model_CNN/plots/ISM/nearbyInPromoter/all_subtypes/nominal_genes.txt")
