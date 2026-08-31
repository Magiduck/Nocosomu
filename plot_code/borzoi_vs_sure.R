library(tidyr)
library(data.table)
library(magrittr)
library(glue)
library(ggplot2)

gene_annotation = fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/shared/ensembl/downloads/genes_Ensembl94_protein_coding.txt")

analysis_id = "nearbyInPromoter"
analysis_id = "alternativeAlleles"

switch = "_abs"
switch = ""
project_folder = "/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/"
bed_file = "1200_upstream_300_downstream_regions_of_protein_coding_genes"
dataset_file = "all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS"
model_file = "SuRE_ResNet_Attention_SuRE4n_tss_selection_m300_p100_stranded_EnhA_intersection_intersection_K562_TSS_EnhA_strand_False_trial_20_model_CNN"
mode_id = glue("ISM{switch}")

analysis_path = glue("/{project_folder}/{bed_file}/{dataset_file}/{model_file}/plots/")

all_sure_tests <- unique(fread(glue("{analysis_path}/{mode_id}/{analysis_id}/all_subtypes/{mode_id}_{analysis_id}_all_subtypes_vs_frequency.tsv")))
all_sure_mutations <- unique(fread(cmd = glue("cat {analysis_path}/{mode_id}/{analysis_id}/all_subtypes/plotData*.tsv | grep -v dummy"), header = F))

project_folder = "/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/"
bed_file = "1200_upstream_300_downstream_regions_of_protein_coding_genes"
dataset_file = "all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS"
mode_id = glue("mutatedBasePairs{switch}")
model_file = "borzoi_gtexBlood"

analysis_path = glue("/{project_folder}/{bed_file}/{dataset_file}/{model_file}/plots/")

all_borzoi_tests <- unique(fread(glue("cat {analysis_path}/{mode_id}/{analysis_id}/all_subtypes/{mode_id}_{analysis_id}_all_subtypes_vs_frequency.tsv")))
all_borzoi_mutations <- unique(fread(cmd = glue("cat {analysis_path}/{mode_id}/{analysis_id}/all_subtypes/plotData*.tsv | grep -v dummy"), header = F))

theshold_show = 20
merged_tests <- all_sure_tests[all_borzoi_tests, on="associated_gene==associated_gene", nomatch = 0]
merged_tests[gene_annotation, on = "associated_gene==V1", symbol := V7]
merged_tests[-log10(i.w_pvalue) > theshold_show | -log10(w_pvalue) > theshold_show, plot_label := symbol]

merged_tests[,cor(-log10(i.w_pvalue),-log10(w_pvalue), method = "spearman")]
ggplot(merged_tests, aes(-log10(w_pvalue), -log10(i.w_pvalue))) + 
  geom_point() +
  geom_text(aes(label = plot_label), hjust=0, vjust = 0)+
  coord_fixed() +
  theme_bw() +
  xlim(0,32) +
  xlab("SuRE ResNet K562") +
  ylab("Borzoi GTEX:blood") +
  ggtitle(glue("{analysis_id} (All subtypes)"))


#how are mutations overall recovered?
all_merged_mutations = all_sure_mutations[all_borzoi_mutations, on = "V1==V1", nomatch = 0]
ggplot(all_merged_mutations, aes(V2, i.V2)) + geom_point(alpha = 0.1) +
  xlab("SuRE ResNet") +
  ylab("Borzoi GTEX:blood") +
  theme_bw() +
  facet_wrap(~V3)

all_merged_mutations[V3=="observed_mutation",cor(V2, i.V2, method = "spearman")]
all_merged_mutations[V3=="comparison_group",cor(V2, i.V2, method = "spearman")]

all_mutations <- fread(cmd = glue("cat {analysis_path}/{mode_id}/{analysis_id}/all_subtypes/plotData*.tsv | grep -v dummy"), header = F)
%>%
  separate(., V1, c("#CHROM", "STA", "POS", "REF", "ALT", "gene"), sep = "_", remove = F)

#get parquet file
parquet_file_mutations = glue("/{project_folder}/{bed_file}/{bed_file}.parquet")


all_mutations <- unique(fread(cmd = glue("cat {analysis_path}/alternativeAlleles/all_subtypes/genes/data/plotData*.tsv {analysis_path}/nearbyInPromoter/all_subtypes/genes/data/plotData*.tsv | grep -v dummy | cut -f1"), header = F)) %>%
  separate(.,
         V1,
         c("#CHROM", "STA", "POS", "REF", "ALT", "gene"), sep = "_", remove = F)

setnames(all_mutations, "V1", "ID")

fwrite(all_mutations[,.(`#CHROM`, POS, ID, REF, ALT)],
       glue("{analysis_path}/significant_nearbyInPromoter_alternativeAlleles_background_mutations.tsv"),
       sep = "\t")


all_mutations <- unique(fread(cmd = glue("cat {analysis_path}/*/all_subtypes/genes/data/plotData*.tsv | grep observ | cut -f1"), header = F)) %>%
  separate(.,
           V1,
           c("#CHROM", "STA", "POS", "REF", "ALT"), sep = "_", remove = F)

setnames(all_mutations, "V1", "ID")

fwrite(all_mutations[,.(`#CHROM`, POS, ID, REF, ALT)],
       glue("{analysis_path}/significant_genes_observed_analyzed_mutations.tsv"),
       sep = "\t")

#all SNP to gene
all_mutations <- fread(cmd = glue("cat {analysis_path}/*/all_subtypes/genes/data/plotData*.tsv | grep -v dummy"), header = F)
setnames(all_mutations, "V1", "ID")
setnames(all_mutations, "V2", "SuRE-ResNet")
setnames(all_mutations, "V5", "gene")

fwrite(unique(all_mutations[,.(`ID`, gene)]),
       glue("{analysis_path}/significant_genes_to_mutations.tsv"),
       sep = "\t")


borzoi_predictions <- fread(glue("/{analysis_path}/significant_genes_observed_analyzed_mutations-borzoi_RNA-blood.tsv"))
borzoi_background_predictions <- fread("/groups/umcg-fg/tmp01/projects/non-coding-somatic/umcg-curzua/pipeline_output/nocosomu_iDriver/6000_upstream_300_downstream_regions_of_protein_coding_genes/all_simple_somatic_mutations_ICGC_release28_chrsorted_WGS/SuRE_ResNet_Attention_SuRE4n_tss_selection_m300_p100_stranded_EnhA_intersection_intersection_K562_TSS_EnhA_strand_False_trial_20_model_CNN/plots/ISM_abs/significant_nearbyInPromoter_alternativeAlleles_background_mutations-borzoi_RNA-all.tsv")
all_mutations[borzoi_predictions, on=.(ID==variant_id), borzoi_logSED := logSED]
all_mutations[borzoi_background_predictions, on=.(ID==variant_id), borzoi_logSED := logSED]


plot_mutations <- all_mutations[V3 == "comparison_group" & !is.na(borzoi_logSED)]

corr_pearson_df <- plot_mutations[!is.na(borzoi_logSED),.SD[, cor(`SuRE-ResNet`, borzoi_logSED, method = "pearson")  %>% round(2)], by=.(gene)] %>% .[order(-abs(V1))]
corr_spearman <- plot_mutations[!is.na(borzoi_logSED), cor(`SuRE-ResNet`, borzoi_logSED, method = "spearman")]  %>% round(2)
ggplot(plot_mutations[!is.na(borzoi_logSED)], aes(`SuRE-ResNet`, borzoi_logSED)) +
  geom_point() +
  ggtitle(glue("Spearman rho: {corr_spearman}")) +
  theme_bw()
ggplot(plot_mutations[!is.na(borzoi_logSED) & gene %in% corr_pearson_df[1:10, gene]], aes(`SuRE-ResNet`, borzoi_logSED)) +
  geom_point() +
  ggtitle(glue("Top 10 genes by Pearson")) +
  theme_bw() +
  facet_wrap(~gene)

corr_spearman <- plot_mutations[!is.na(borzoi_logSED) & gene %in% corr_pearson_df[1:10, gene], cor(`SuRE-ResNet`, borzoi_logSED, method = "spearman")]  %>% round(2)
ggplot(plot_mutations[!is.na(borzoi_logSED) & gene %in% corr_pearson_df[1:10, gene]], aes(`SuRE-ResNet`, borzoi_logSED)) +
  geom_point() +
  ggtitle(glue("Top 10 genes by Pearson")) +
  facet_wrap(~gene) +
  theme_bw() + ylim(-0.2,0.2)

ggplot(plot_mutations[!is.na(borzoi_logSED) & gene %in% corr_pearson_df[55:65, gene]], aes(`SuRE-ResNet`, borzoi_logSED)) +
  geom_point() +
  ggtitle(glue("Top 10 genes by Pearson")) +
  theme_bw()

ggplot(all_mutations[!is.na(borzoi_logSED) & gene == "ENSG00000162062"], aes(borzoi_logSED)) +
  geom_histogram(aes(col=V3)) +
  ggtitle(glue("TERT")) +
  theme_bw()

https://geneticsumcg.slack.com/archives/D037H1GC2RW/p1705053849766139
